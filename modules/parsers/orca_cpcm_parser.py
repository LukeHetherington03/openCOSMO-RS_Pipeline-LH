class OrcaCpcmParser:
    def __init__(self, filepath):
        self.filepath = filepath

    def parse(self):
        data = {}
        lines = open(self.filepath).read().splitlines()

        # Basic fields
        data["n_atoms"] = int(lines[0].split()[0])
        data["n_segments"] = int(lines[1].split()[0])
        data["surface_type"] = int(lines[2].split()[0])
        data["epsilon_function_type"] = int(lines[3].split()[0])
        data["print_level"] = int(lines[4].split()[0])
        data["feps_x_flag"] = int(lines[5].split()[0])
        data["lebedev_points"] = int(lines[7].split()[0])
        data["isodensity_scheme"] = int(lines[9].split()[0])

        # Thresholds
        data["threshold_h"] = float(lines[11].split()[0])
        data["threshold_non_h"] = float(lines[12].split()[0])

        # FEps X parameter
        data["feps_x_parameter"] = float(lines[14].split()[0])

        # Cutoffs
        data["cutoff_area"] = float(lines[16].split()[0])
        data["cutoff_switch"] = float(lines[17].split()[0])

        # Volume + Area
        data["volume"] = float(lines[19].split()[0])
        data["area"] = float(lines[20].split()[0])

        # Energies
        data["cpcm_dielectric_energy"] = float(lines[22].split()[0])
        data["one_electron_energy"] = float(lines[23].split()[0])

        return data
