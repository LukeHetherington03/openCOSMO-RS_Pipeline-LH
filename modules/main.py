from .work_flow_manager import WorkflowManager

def main():
    wm = WorkflowManager()

    # Resume from pruning stage only
    wm.run_pipeline(
        pruning_args={
            "method_subdir": "acrylates/rdkit",
            "method": "topN",
            "top_n": 5
        },
        optimisation_args={
            "engine": "orca",
            "functional": "B3LYP",
            "basis": "def2-SVP"
        },
    )


# ------------------ Run Pipeline ------------------

if __name__ == "__main__":
    main()
