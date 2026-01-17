import os
import json
import subprocess
import shutil
from pathlib import Path

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter


class SolubilityStage(BaseStage):
    """
    Solubility stage using openCOSMO-RS:
      - Reads orcacosmo_summary.json
      - Creates mixture input files for openCOSMO-RS
      - Runs solubility calculations
      - Parses results into solubility_summary.json
    """

    def execute(self):
        self.log_header("Starting Solubility Stage")

        # Get openCOSMO-RS binary
        self.opencosmo_binary = self.parameters.get("opencosmo_binary", "openCOSMORS")

        # Solubility calculation parameters
        self.temperature = self.parameters.get("temperature", 298.15)
        self.solvent_name = self.parameters.get("solvent_name", "water")
        self.solvent_cosmo_dir = self.parameters.get("solvent_cosmo_dir", "openCOSMO-RS_Pipeline-LH/CONSTANT_FILES")
        self.solvent_smiles = self.parameters.get("solvent_smiles", "O")
        
        # Melting point parameters for solubility iteration
        self.melting_temp = self.parameters.get("melting_temp", "N/A")
        self.Gfus = self.parameters.get("Gfus", "N/A")
        self.Hfus = self.parameters.get("Hfus", "N/A")
        self.SORcf = self.parameters.get("SORcf", 1.0)
        
        # Calculation type
        self.calculations = self.parameters.get("calculations", "all")  # 'all', 'pure_only', 'mixed_only'

        if not self.solvent_cosmo_dir:
            self.fail("solvent_cosmo_dir is required for solubility calculations")

        if not os.path.exists(self.solvent_cosmo_dir):
            self.fail(f"Solvent COSMO directory not found: {self.solvent_cosmo_dir}")

        # Prepare directories
        self._prepare_directories()

        # Load orcacosmo summary
        self._load_orcacosmo_summary()

        # Discover solutes from summary
        solute_ids = [entry["lookup_id"] for entry in self.orcacosmo_entries]
        self.set_items(solute_ids)

        # Track successful calculations
        self.successful_calcs = []

        # Process each solute
        for lookup_id in list(self.job.pending_items):
            try:
                self._process_solute(lookup_id)
                self.log(f"Completed {lookup_id}")
                self.update_progress(lookup_id)
            except Exception as e:
                self.log(f"[ERROR] {lookup_id}: {e}")
                self.update_progress(lookup_id, success=False)

        # Write summary
        self._write_summary()

        self.job.mark_complete()
        self.log_header("Solubility Stage Complete")

    # ------------------------------------------------------------
    # Directory setup
    # ------------------------------------------------------------
    def _prepare_directories(self):
        self.workdir = os.path.join(self.outputs_dir, "opencosmo_workdir")
        self.results_dir = os.path.join(self.outputs_dir, "results")
        self.cosmo_copies_dir = os.path.join(self.outputs_dir, "cosmo_files")

        os.makedirs(self.workdir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.cosmo_copies_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Load orcacosmo summary
    # ------------------------------------------------------------
    def _load_orcacosmo_summary(self):
        summary_file = self.parameters.get("summary_file")
        if not summary_file:
            self.fail("Solubility stage requires summary_file")

        if not os.path.exists(summary_file):
            self.fail(f"summary_file does not exist: {summary_file}")

        with open(summary_file) as f:
            self.orcacosmo_entries = json.load(f)

        self.log(f"Loaded {len(self.orcacosmo_entries)} entries from orcacosmo summary")

    # ------------------------------------------------------------
    # Process single solute
    # ------------------------------------------------------------
    def _process_solute(self, lookup_id):
        # Find entry in summary
        entry = next((e for e in self.orcacosmo_entries if e["lookup_id"] == lookup_id), None)
        if not entry:
            raise ValueError(f"No entry found for {lookup_id}")

        orcacosmo_path = entry["orcacosmo_path"]
        if not os.path.exists(orcacosmo_path):
            raise FileNotFoundError(f"ORCACOSMO file missing: {orcacosmo_path}")

        # Copy .orcacosmo file to a dedicated directory structure
        solute_cosmo_dir = os.path.join(self.cosmo_copies_dir, lookup_id)
        os.makedirs(solute_cosmo_dir, exist_ok=True)
        
        dest_cosmo = os.path.join(solute_cosmo_dir, f"{lookup_id}.orcacosmo")
        shutil.copy(orcacosmo_path, dest_cosmo)

        # Get solute SMILES if available (optional)
        solute_smiles = entry.get("smiles", "")

        # Create openCOSMO-RS input file
        input_file = os.path.join(self.workdir, f"{lookup_id}_input.txt")
        self._write_opencosmo_input(
            input_file,
            lookup_id,
            solute_cosmo_dir,
            solute_smiles
        )

        # Run openCOSMO-RS
        output_file = os.path.join(self.workdir, f"{lookup_id}_output.txt")
        self._run_opencosmo(input_file, output_file)

        # Parse results
        results = self._parse_output(output_file, lookup_id)

        # Save results JSON
        results_json = os.path.join(self.results_dir, f"{lookup_id}_solubility.json")
        with AtomicWriter(results_json) as f:
            json.dump(results, f, indent=2)

        # Track success
        self.successful_calcs.append({
            "lookup_id": lookup_id,
            "results_path": results_json,
            "orcacosmo_path": orcacosmo_path,
            **results
        })

    # ------------------------------------------------------------
    # Write openCOSMO-RS input file
    # ------------------------------------------------------------
    def _write_opencosmo_input(self, path, solute_name, solute_cosmo_dir, solute_smiles):
        """
        Write the openCOSMO-RS input format.
        """
        # Determine saturation parameter
        if self.melting_temp != "N/A":
            saturation_line = f"saturation {solute_name}"
            if solute_smiles:
                saturation_line += f" {solute_smiles}"
        else:
            saturation_line = "saturation no"

        # Count number of .orcacosmo files in solute directory
        n_solute_confs = len(list(Path(solute_cosmo_dir).glob("*.orcacosmo")))
        
        # Count solvent conformers
        n_solvent_confs = len(list(Path(self.solvent_cosmo_dir).glob("*.orcacosmo")))

        # Generate multiplicity strings (all 1s)
        solute_multiplicities = " ".join(["1"] * n_solute_confs)
        solvent_multiplicities = " ".join(["1"] * n_solvent_confs)

        with open(path, "w") as f:
            f.write(f"calculations {self.calculations}\n")
            f.write(f"temperature {self.temperature}\n")
            f.write(f"{saturation_line}\n")
            f.write(f"meltingtemp {self.melting_temp}\n")
            f.write(f"Gfus {self.Gfus}\n")
            f.write(f"Hfus {self.Hfus}\n")
            f.write(f"SORcf {self.SORcf}\n")
            f.write("# name # mol fraction (of molecule) # path_to_dir # nconf # multiplicities\n")
            
            # Solute line
            f.write(f"{solute_name}\t0.0\t{solute_cosmo_dir}\t{n_solute_confs}\t{solute_multiplicities}\n")
            
            # Solvent line
            f.write(f"{self.solvent_name}\t1.0\t{self.solvent_cosmo_dir}\t{n_solvent_confs}\t{solvent_multiplicities}\n")

    # ------------------------------------------------------------
    # Run openCOSMO-RS
    # ------------------------------------------------------------
    def _run_opencosmo(self, input_file, output_file):
        """
        Execute the openCOSMO-RS binary.
        """
        with open(output_file, "w") as f:
            result = subprocess.run(
                [self.opencosmo_binary, input_file],
                stdout=f,
                stderr=subprocess.STDOUT,
                check=False  # Don't raise on non-zero exit
            )

        if result.returncode != 0:
            # Log warning but don't fail - output file may still have useful info
            self.log(f"[WARNING] openCOSMO-RS returned code {result.returncode}")

    # ------------------------------------------------------------
    # Parse openCOSMO-RS output
    # ------------------------------------------------------------
    def _parse_output(self, output_file, lookup_id):
        """
        Parse openCOSMO-RS output file.
        This is a basic parser - adjust based on actual output format.
        """
        results = {
            "lookup_id": lookup_id,
            "temperature": self.temperature,
            "solvent": self.solvent_name,
            "raw_output": None,
            "activity_coefficient": None,
            "solubility_mol_frac": None,
            "solubility_g_L": None,
        }

        if not os.path.exists(output_file):
            return results

        with open(output_file) as f:
            output_text = f.read()
            results["raw_output"] = output_text

        # Parse key values from output
        # TODO: Adjust these parsers based on actual openCOSMO-RS output format
        for line in output_text.split('\n'):
            line = line.strip()
            
            # Example parsing - adjust to match actual output
            if "activity coefficient" in line.lower():
                try:
                    results["activity_coefficient"] = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
            
            if "solubility" in line.lower() and "mol/mol" in line.lower():
                try:
                    results["solubility_mol_frac"] = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
            
            if "solubility" in line.lower() and "g/L" in line.lower():
                try:
                    results["solubility_g_L"] = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass

        return results

    # ------------------------------------------------------------
    # Write summary
    # ------------------------------------------------------------
    def _write_summary(self):
        """
        Write solubility_summary.json for potential downstream stages.
        """
        summary_path = os.path.join(self.outputs_dir, "solubility_summary.json")
        
        with AtomicWriter(summary_path) as f:
            json.dump(self.successful_calcs, f, indent=2)
        
        self.log(f"Wrote summary with {len(self.successful_calcs)} entries: {summary_path}")