#!/usr/bin/env python
"""
05_run_md.py — Short explicit-solvent MD of one CDK2-ligand complex with OpenMM (GPU).

Reads the Amber topology built by 04_setup_complex.py and runs:
  minimize -> NVT heat to 300 K (heavy-atom restrained) -> NPT equilibration (restraints
  released) -> production (NPT, 300 K).  The production trajectory is saved as DCD and as an
  Amber NetCDF (for MM-GBSA), and the ligand heavy-atom RMSD-to-start is written out as the
  pose-stability metric (docked poses that drift away are flagged).

Outputs (under <workdir>/md/<cid>/):
  traj.dcd, traj.nc        production trajectory (gitignored)
  equilibrated.pdb         topology/first frame for analysis
  rmsd.csv                 time(ns), ligand_rmsd(A), protein_bb_rmsd(A)
  md_summary.json          mean/final ligand RMSD, stable flag, ns, n_frames
"""
import argparse
import json
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit


def log(msg):
    print(f"[md] {msg}", flush=True)


def restraint_force(prmtop, positions, k):
    force = mm.CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    for atom in prmtop.topology.atoms():
        if atom.residue.name not in ("HOH", "WAT", "Na+", "Cl-", "K+") and atom.element != app.element.hydrogen:
            force.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometer))
    return force


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--cid", required=True)
    ap.add_argument("--ns", type=float, default=4.0, help="production length (ns)")
    ap.add_argument("--equil-ps", type=float, default=300.0)
    ap.add_argument("--report-ps", type=float, default=10.0)
    ap.add_argument("--temp", type=float, default=300.0)
    args = ap.parse_args()

    wd = Path(args.workdir).resolve()
    d = wd / "md" / args.cid
    prmtop = app.AmberPrmtopFile(str(d / "complex_solv.prmtop"))
    inpcrd = app.AmberInpcrdFile(str(d / "complex_solv.inpcrd"))

    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                                 constraints=app.HBonds)
    dt = 0.002 * unit.picoseconds
    integrator = mm.LangevinMiddleIntegrator(args.temp * unit.kelvin, 1.0 / unit.picosecond, dt)
    platform = mm.Platform.getPlatformByName("CUDA")
    props = {"Precision": "mixed"}

    # restrained equilibration
    restr = restraint_force(prmtop, inpcrd.positions, 5.0 * unit.kilocalories_per_mole / unit.angstrom**2)
    restr_idx = system.addForce(restr)
    sim = app.Simulation(prmtop.topology, system, integrator, platform, props)
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    log("minimizing ...")
    sim.minimizeEnergy(maxIterations=5000)

    log(f"NVT heat + restrained equilibration ({args.equil_ps} ps) ...")
    sim.context.setVelocitiesToTemperature(args.temp * unit.kelvin)
    sim.step(int((args.equil_ps * 0.3) / 0.002))  # NVT, restrained
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, args.temp * unit.kelvin))
    sim.context.reinitialize(preserveState=True)
    sim.step(int((args.equil_ps * 0.3) / 0.002))  # NPT, still restrained
    # release restraints gradually
    for k in (2.0, 0.5, 0.0):
        sim.context.setParameter("k", k * 4.184 * 100)  # kcal/mol/A^2 -> kJ/mol/nm^2
        sim.step(int((args.equil_ps * 0.4 / 3) / 0.002))

    # production
    n_steps = int(args.ns * 1000 / 0.002)
    report = int(args.report_ps / 0.002)
    sim.reporters.append(app.DCDReporter(str(d / "traj.dcd"), report))
    sim.reporters.append(app.StateDataReporter(
        str(d / "md.log"), report * 5, step=True, time=True, potentialEnergy=True,
        temperature=True, density=True, speed=True))
    with open(d / "equilibrated.pdb", "w") as fh:
        app.PDBFile.writeFile(prmtop.topology,
                              sim.context.getState(getPositions=True).getPositions(), fh)
    log(f"production {args.ns} ns ({n_steps} steps), reporting every {args.report_ps} ps ...")
    sim.step(n_steps)
    log("production done; analyzing RMSD")

    analyze(d, args.report_ps)


def analyze(d, report_ps):
    import mdtraj as md
    traj = md.load(str(d / "traj.dcd"), top=str(d / "equilibrated.pdb"))
    traj.save_netcdf(str(d / "traj.nc"))  # for MM-GBSA
    prot_bb = traj.topology.select("protein and backbone")
    lig = traj.topology.select("resname LIG and not element H")
    ref = traj[0]
    traj.superpose(ref, atom_indices=prot_bb)
    lig_rmsd = np.sqrt(np.mean(np.sum((traj.xyz[:, lig] - ref.xyz[0, lig])**2, axis=2), axis=1)) * 10.0
    bb_rmsd = md.rmsd(traj, ref, atom_indices=prot_bb) * 10.0
    times = np.arange(len(traj)) * report_ps / 1000.0
    np.savetxt(d / "rmsd.csv", np.column_stack([times, lig_rmsd, bb_rmsd]),
               delimiter=",", header="time_ns,ligand_rmsd_A,protein_bb_rmsd_A", comments="")
    summary = {
        "cid": d.name, "n_frames": int(len(traj)), "ns": float(times[-1]) if len(times) else 0.0,
        "ligand_rmsd_mean": round(float(lig_rmsd.mean()), 2),
        "ligand_rmsd_final": round(float(lig_rmsd[-1]), 2),
        "ligand_rmsd_max": round(float(lig_rmsd.max()), 2),
        "protein_bb_rmsd_mean": round(float(bb_rmsd.mean()), 2),
        "stable": bool(lig_rmsd.mean() < 3.0 and lig_rmsd[-1] < 4.0),
    }
    (d / "md_summary.json").write_text(json.dumps(summary, indent=2))
    log(f"summary: {summary}")


if __name__ == "__main__":
    main()
