from work_flow_manager import WorkflowManager

def main():
    manager = WorkflowManager()

    # ------------------ Step 1: Clean raw dataset ------------------
    # inp_file = step_1(manager)

    # ------------------ Step 2: Generate conformers ------------------
    # inp_file_2 = "acrylates/acrylate_initiator_set_clean.inp"
    # xyz_dir = step_2(manager, inp_file_2)

    # ------------------ Step 3: Prune conformers ------------------
    #pruned_lookup = step_3(manager, "acrylates/rdkit", method="topN", top_n=3)
    #print(f"Pruned lookup file: {pruned_lookup}")

    # ------------------ Step 4: Optimise conformers ------------------
    #optimisation_summary = step_4(manager, "acrylates/rdkit/topN_3", engine="gxtb", level="GFN2-xTB")
    #print(f"Optimisation summary written to: {optimisation_summary}")

    # ------------------ Step 5: Generate COSMO files ------------------
    mixture_inputs = step_5(manager, "acrylates/rdkit/topN_5/gxtb_GFN2-xTB")
    #print(f"Mixture input file ready at: {mixture_inputs}")


# ------------------ Step Definitions ------------------

def step_1(manager, raw_filename="acrylate_initiator_set.csv"):
    """Step 1: Clean raw dataset."""
    inp_file = manager.run_data_cleaning(raw_filename)
    print(f"Conformer input ready at {inp_file}")
    return inp_file


def step_2(manager, inp_file_relpath):
    """Step 2: Generate conformers with RDKit defaults."""
    generation_statement = {
        "engine": "rdkit",
        "n_conformers": 250,
        "seed": 42
    }
    xyz_dir = manager.run_conformer_generation(inp_file_relpath, generation_statement)
    print(f"Conformers written to: {xyz_dir}")
    return xyz_dir


def step_3(manager, method_subdir, method="energy_window", **params):
    """
    Step 3: Prune conformers.
    Uses WorkflowManager to dispatch pruning strategies.
    Output folder follows: 4_pruned_conformers/<dataset>/<engine>/<method>_<param>
    """
    return manager.run_conformer_pruning(method_subdir, method=method, **params)


def step_4(manager, pruned_subdir, engine="gxtb", **params):
    """
    Step 4: Optimise conformers.
    Uses WorkflowManager to run geometry optimisation with xTB or ORCA.
    Output folder follows:
        5_conformer_xyz_optimised/<dataset>/<engine>/<method>_<param>/<engine_engineparam>
    """
    return manager.run_conformer_optimisation(pruned_subdir, engine=engine, **params)

def step_5(manager, optimisation_subdir,
           method="B3LYP", basis="def2-SVP",
           solvent="Water", charge=0, multiplicity=1):
    print("=== Step 5: ORCA COSMO runs ===")
    results = manager.run_orca_cosmo_step(
        optimisation_subdir,
        method=method,
        basis=basis,
        solvent=solvent,
        charge=charge,
        multiplicity=multiplicity
    )
    for inchi, info in results.items():
        print(f"{inchi}: {info['n_confs']} conformers → {info['dir']}")
    return results


# ------------------ Run Pipeline ------------------

if __name__ == "__main__":
    main()
