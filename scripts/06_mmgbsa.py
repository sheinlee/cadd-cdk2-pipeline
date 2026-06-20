#!/usr/bin/env python
"""
06_mmgbsa.py — MM-GBSA binding free-energy rescoring of one MD trajectory (AmberTools MMPBSA.py).

Single-trajectory MM-GBSA: the solvated production trajectory is stripped to complex/receptor/
ligand and the GB binding free energy is averaged over frames. Optional per-residue energy
decomposition (over the residues lining the pocket) identifies the binding hotspots.

If the decomposition run fails for any reason, MM-GBSA is retried without it so the core
dG_bind is still produced.

Outputs (under <workdir>/md/<cid>/):
  FINAL_RESULTS_MMPBSA.dat   full MMPBSA report
  FINAL_DECOMP_MMPBSA.dat    per-residue decomposition (if available)
  mmgbsa.json                parsed dG_bind (mean/std/sem) + top contributing residues
"""
import argparse
import json
import re
import subprocess
from pathlib import Path


def log(msg):
    print(f"[mmgbsa] {msg}", flush=True)


def pocket_residues(d, cutoff_nm=0.6):
    """1-based complex-prmtop residue numbers within `cutoff` of the ligand (+ the ligand)."""
    import mdtraj as md
    t = md.load(str(d / "complex_dry.inpcrd"), top=str(d / "complex_dry.prmtop"))
    lig = t.topology.select("resname LIG")
    prot = t.topology.select("protein")
    near = md.compute_neighbors(t, cutoff_nm, lig, haystack_indices=prot)[0]
    res = {t.topology.atom(i).residue.index + 1 for i in near}
    res.add(t.topology.atom(lig[0]).residue.index + 1)
    return ",".join(str(r) for r in sorted(res))


def write_input(path, interval, print_res):
    txt = f"""MM-GBSA for CDK2-ligand complex
&general
  startframe=1, interval={interval}, verbose=2, keep_files=0,
/
&gb
  igb=2, saltcon=0.150,
/
"""
    if print_res:
        txt += f"""&decomp
  idecomp=1, dec_verbose=0, print_res="{print_res}",
/
"""
    path.write_text(txt)


def run_mmpbsa(d, interval, print_res):
    write_input(d / "mmpbsa.in", interval, print_res)
    cmd = ["MMPBSA.py", "-O", "-i", "mmpbsa.in", "-o", "FINAL_RESULTS_MMPBSA.dat"]
    if print_res:
        cmd += ["-do", "FINAL_DECOMP_MMPBSA.dat"]
    cmd += ["-sp", "complex_solv.prmtop", "-cp", "complex_dry.prmtop",
            "-rp", "receptor_dry.prmtop", "-lp", "ligand_dry.prmtop", "-y", "traj.nc"]
    log("$ " + " ".join(cmd))
    return subprocess.run(cmd, cwd=d, capture_output=True, text=True)


def parse_results(dat):
    mean = std = sem = None
    for line in dat.read_text().splitlines():
        if line.strip().startswith("DELTA TOTAL"):
            nums = re.findall(r"[-+]?\d+\.\d+", line)
            if len(nums) >= 3:
                mean, std, sem = float(nums[0]), float(nums[1]), float(nums[2])
    return mean, std, sem


def parse_decomp(dat, top_n=8):
    if not dat.exists():
        return []
    res, started = [], False
    for line in dat.read_text().splitlines():
        if line.strip().startswith("Resid"):
            started = True
            continue
        if started:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 18:
                try:
                    res.append((parts[0], float(parts[-3])))
                except ValueError:
                    pass
    res.sort(key=lambda x: x[1])
    return [{"residue": r, "kcal_mol": round(v, 2)} for r, v in res[:top_n]]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--cid", required=True)
    ap.add_argument("--interval", type=int, default=4)
    ap.add_argument("--decomp", action="store_true", default=True)
    ap.add_argument("--no-decomp", dest="decomp", action="store_false")
    args = ap.parse_args()

    wd = Path(args.workdir).resolve()
    d = wd / "md" / args.cid

    print_res = None
    if args.decomp:
        try:
            print_res = pocket_residues(d)
            log(f"pocket residues for decomposition: {print_res}")
        except Exception as exc:
            log(f"could not determine pocket residues ({exc.__class__.__name__}); skipping decomp")

    r = run_mmpbsa(d, args.interval, print_res)
    decomp_done = print_res is not None
    if r.returncode != 0 and print_res is not None:
        log("decomposition run failed; retrying MM-GBSA without per-residue decomposition")
        r = run_mmpbsa(d, args.interval, None)
        decomp_done = False
    if r.returncode != 0:
        raise RuntimeError(f"MMPBSA.py failed (rc={r.returncode})\n{r.stdout[-2000:]}\n{r.stderr[-1500:]}")

    mean, std, sem = parse_results(d / "FINAL_RESULTS_MMPBSA.dat")
    out = {
        "cid": args.cid,
        "dG_bind_kcal_mol": None if mean is None else round(mean, 2),
        "std": None if std is None else round(std, 2),
        "sem": None if sem is None else round(sem, 2),
        "top_residues": parse_decomp(d / "FINAL_DECOMP_MMPBSA.dat") if decomp_done else [],
    }
    (d / "mmgbsa.json").write_text(json.dumps(out, indent=2))
    log(f"result: {out}")


if __name__ == "__main__":
    main()
