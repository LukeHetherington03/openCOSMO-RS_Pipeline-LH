from .orca_log_parser_rev import OrcaLogParserRev
from .orca_cpcm_parser_rev import OrcaCpcmParserRev
from .orca_cpcm_corr_parser_rev import OrcaCpcmCorrParserRev

from modules.utils.molecule_utils import MoleculeUtils


import json

class OrcaCosmoOrchestrator:
    def __init__(self, bundle):
        self.meta = bundle["meta"]
        self.paths = bundle["paths"]

        # Load XYZ
        self.xyz_atoms = self._load_xyz(self.paths["xyz"])

        # Build adjacency
        self.adjacency = MoleculeUtils.adjacency_from_xyz(self.xyz_atoms)

        # ------------------------------------------------------------
        # FIX: Load JSON from disk before passing to reverse parsers
        # ------------------------------------------------------------
        with open(self.paths["log"]) as f:
            log_json = json.load(f)

        with open(self.paths["cpcm"]) as f:
            cpcm_json = json.load(f)

        with open(self.paths["cpcm_corr"]) as f:
            cpcm_corr_json = json.load(f)

        # Instantiate reverse parsers with JSON, not paths
        self.log_parser = OrcaLogParserRev(log_json)
        self.cpcm_parser = OrcaCpcmParserRev(cpcm_json)
        self.cpcm_corr_parser = OrcaCpcmCorrParserRev(cpcm_corr_json)

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def reconstruct(self):
        out = []

        # Canonical ORCA-COSMO block order
        out.append(self._header_block())
        out.append(self.log_parser.reconstruct_energy())
        out.append(self.log_parser.reconstruct_dipole())
        out.append(self._xyz_block())

        out.append("\n##################################################")
        out.append("#COSMO")
        out.append(self.cpcm_parser.reconstruct())

        out.append(self.cpcm_corr_parser.reconstruct())
        out.append(self._adjacency_block())
        out.append(self.log_parser.reconstruct_polarizabilities())
        # LOG blocks in correct order


        return "\n".join(out)

    # ------------------------------------------------------------
    # Header
    # ------------------------------------------------------------
    def _header_block(self):
        m = self.meta
        return (
            f"{m['lookup_id']} : DFT_CPCM_BP86_{m['method_used']}\n"
        )

    # ------------------------------------------------------------
    # XYZ loader
    # ------------------------------------------------------------
    def _load_xyz(self, path):
        atoms = []
        with open(path) as f:
            lines = f.readlines()

        start = 2 if len(lines) > 2 and lines[0].strip().isdigit() else 0

        for line in lines[start:]:
            parts = line.split()
            if len(parts) == 4:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))

        return atoms

    # ------------------------------------------------------------
    # XYZ block
    # ------------------------------------------------------------
    def _xyz_block(self):
        out = []
        atoms = self.xyz_atoms

        out.append("##################################################")
        out.append("#XYZ_FILE")
        out.append(str(len(atoms)))

        # Energy comes from the log parser's parsed JSON
        energy = self.log_parser.data.get("final_single_point_energy")

        out.append(
            f"Coordinates from ORCA-job {self.meta['method_used']} "
            f"E {energy}"
        )

        for el, x, y, z in atoms:
            out.append(f"  {el:<2}  {x:20.14f}  {y:20.14f}  {z:20.14f}")

        return "\n".join(out)

    # ------------------------------------------------------------
    # Adjacency block
    # ------------------------------------------------------------
    def _adjacency_block(self):
        out = []
        out.append("##################################################")
        out.append("#ADJACENCY_MATRIX")

        for row in self.adjacency:
            out.append("".join(f"{int(v):4d}" for v in row))

        return "\n".join(out)
