#!/usr/bin/env python3

import os
import shutil
import subprocess
import json

# -------------------------------------------------------------------
# Paths (adjust if needed)
# -------------------------------------------------------------------
XYZ = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/requests/R-20Jan26-1025-05841-acr_t2_test/jobs/J-20Jan26-1028-23611-orcacosmo/inputs/BAPJBEWLBFYGME-UHFFFAOYSA-N_conf008.xyz"
ORCA = "/home/lunet/cglh4/software/orca_6_0_0/orca"
OUTROOT = "orca_variants"

CPCM_JSON = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/CONSTANT_FILES/chemistry/cpcm_radii.json"

# -------------------------------------------------------------------
# Load full CPCM radii + cut_area from JSON
# -------------------------------------------------------------------
with open(CPCM_JSON) as f:
    cpcm_data = json.load(f)

radii = cpcm_data["radii"]
cut_area = cpcm_data["cut_area"]


def build_cpcm_block() -> str:
    lines = ["%cpcm"]
    for Z, r in radii.items():
        lines.append(f"  radius[{Z}]  {r}")
    lines.append(f"  cut_area {cut_area}")
    lines.append("end")
    return "\n".join(lines)


# -------------------------------------------------------------------
# Define different ORCA single-stage COSMO variants
# -------------------------------------------------------------------
variants = [
    {
        "name": "v01_bp86_def2_tzvpd_elprop",
        "header": "! CPCM BP86 def2-TZVPD SP",
        "elprop": True,
        "cpcm": True,
    },
    {
        "name": "v03_bp86_def2_tzvp_elprop",
        "header": "! CPCM BP86 def2-TZVP SP",
        "elprop": True,
        "cpcm": True,
    },

]

# -------------------------------------------------------------------
# Prepare output root
# -------------------------------------------------------------------
os.makedirs(OUTROOT, exist_ok=True)

# -------------------------------------------------------------------
# Generate and run each variant
# -------------------------------------------------------------------
cpcm_block = build_cpcm_block()

for v in variants:
    vdir = os.path.join(OUTROOT, v["name"])
    os.makedirs(vdir, exist_ok=True)

    # Copy XYZ
    xyz_name = os.path.basename(XYZ)
    xyz_local = os.path.join(vdir, xyz_name)
    shutil.copy(XYZ, xyz_local)

    # Write ORCA input
    inp = os.path.join(vdir, "input.inp")
    with open(inp, "w") as f:
        f.write("%MaxCore 2000\n\n")
        f.write('%base "cosmo_test"\n\n')

        if v["cpcm"]:
            f.write(cpcm_block + "\n\n")

        f.write(v["header"] + "\n\n")

        if v["elprop"]:
            f.write("%elprop\n")
            f.write("  Polar 1\n")
            f.write("  Polaratom 1\n")
            f.write("end\n\n")

        f.write(f"* xyzfile 0 1 {xyz_name}\n")

    # Run ORCA 6.0.0
    log = os.path.join(vdir, "output.log")
    with open(log, "w") as f:
        subprocess.run(
            [ORCA, "input.inp"],
            cwd=vdir,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=False,
        )

print("All ORCA variants generated and executed under:", OUTROOT)
