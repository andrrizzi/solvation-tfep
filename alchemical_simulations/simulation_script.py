import BioSimSpace as BSS

def BSS_simulation(smilesA: str, smilesB: str, durations) -> None :    
    # moleculeA and moleculeB topology from SMILES
    moleculeA = BSS.Parameters.gaff(smilesA).getMolecule() 
    moleculeB = BSS.Parameters.gaff(smilesB).getMolecule() 

    # alignment of molecules through RMSD, since methane and iodomethane they are rigid
    mapping = BSS.Align.matchAtoms(moleculeA, moleculeB)
    moleculeA = BSS.Align.rmsdAlign(moleculeA, moleculeB, mapping)

    # generation of merged topology for alchemical calculation
    merged = BSS.Align.merge(moleculeA, moleculeB, mapping)

    # solvation of system in a cubic box with edge 5nm
    solvated = BSS.Solvent.tip3p(molecule = merged, box = 3 * [5 * BSS.Units.Length.nanometer])

    # gromacs minimisation, equilibration and production
    eq_time, prod_time = durations
    min_protocol = BSS.Protocol.FreeEnergyMinimisation(steps = 5000)
    min_process = BSS.MD.run(solvated, min_protocol)
    minimised = min_process.getSystem(block = True)

    eq_protocol = BSS.Protocol.FreeEnergyEquilibration(runtime = eq_time)
    eq_process = BSS.MD.run(minimised, eq_protocol)
    equilibrated = eq_process.getSystem(block = True)

    prod_protocol = BSS.Protocol.FreeEnergyProduction(runtime = prod_time)
    free_gmx = BSS.FreeEnergy.Relative(equilibrated, prod_protocol, engine = "gromacs", work_dir = "freenrg_gmx/free")
    BSS.FreeEnergy.Relative.run(free_gmx)
    vac_gmx = BSS.FreeEnergy.Relative(merged.toSystem(), prod_protocol, engine = "gromacs", work_dir = "freenrg_gmx/vacuum")
    BSS.FreeEnergy.Relative.run(vac_gmx)

    # processes continue in backgroun once launches, remember to let them conclude
    # before trying to analyse results or you'll get mysterious errors
    free_gmx.wait()
    vac_gmx.wait()

    # results analysis
    pmf, overlap = BSS.FreeEnergy.Relative.analyse("./freenrg_gmx")
    is_okay = BSS.FreeEnergy.Relative.checkOverlap(overlap)[0]
    if is_okay :
        print(BSS.FreeEnergy.Relative.difference(pmf))
    else :
        print("Results are not reliable, please extend duration of the simulation.")

if __name__ == '__main__':
    duration = { "equilibration": 40, "production": 100 } # in picoseconds
    mol_A = "C"
    mol_B = "CI"

    duration = (BSS.Types.Time(duration["equilibration"], "ps"), BSS.Types.Time(duration["production"], "ps"))

    BSS_simulation(mol_A, mol_B, duration)
