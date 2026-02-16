PRUNING STAGE — ARGUMENT REFERENCE
==================================

ARGUMENT LIST (QUICK REFERENCE)
-------------------------------
energy_window
energy_window_units
rmsd_threshold
max_energy
percentile
n
n_high
n_start
keep_all
strict


DETAILED ARGUMENT DESCRIPTIONS
==============================

1. ENERGY-BASED ARGUMENTS
-------------------------

energy_window
    Type: float or null
    Default: 5.0
    Meaning: Keeps all conformers within ΔE of the lowest-energy conformer.
    Units: Controlled by energy_window_units.

energy_window_units
    Type: string
    Default: "kcal"
    Allowed: "kcal", "kJ", "hartree"
    Meaning: Determines the units used for energy_window.

max_energy
    Type: float or null
    Default: null
    Meaning: Keeps only conformers with absolute energy ≤ this value (in Hartree).

percentile
    Type: float or null
    Default: null
    Meaning: Keeps conformers whose energy is below the given percentile of the
             molecule’s energy distribution.


2. GEOMETRY-BASED ARGUMENTS
---------------------------

rmsd_threshold
    Type: float or null
    Default: 0.5
    Meaning: RMSD-based diversity pruning. Conformers within this RMSD of a
             lower-energy conformer are removed.
    Units: Ångström.


3. COUNT-BASED ARGUMENTS
------------------------

n
    Type: integer or null
    Default: null
    Meaning: Keep the lowest-energy N conformers.

n_high
    Type: integer or null
    Default: null
    Meaning: Keep the highest-energy N conformers.

n_start
    Type: integer or null
    Default: null
    Meaning: Starting index for slice pruning. Used together with n to select a
             contiguous block of conformers.


4. BEHAVIOURAL ARGUMENTS
------------------------

keep_all
    Type: boolean
    Default: false
    Meaning: If true, bypasses all pruning and keeps every conformer.

strict
    Type: boolean
    Default: false
    Meaning: If true, the stage fails when all conformers for a molecule are
             removed.


SUMMARY TABLE
=============

Argument              | Type        | Default | Description
--------------------- | ----------- | ------- | -----------------------------------------------
energy_window         | float/null  | 5.0     | ΔE window relative to minimum energy
energy_window_units   | string      | "kcal"  | Units for energy window
rmsd_threshold        | float/null  | 0.5     | RMSD diversity pruning threshold
max_energy            | float/null  | null    | Absolute energy cutoff
percentile            | float/null  | null    | Percentile-based pruning
n                     | int/null    | null    | Keep lowest N conformers
n_high                | int/null    | null    | Keep highest N conformers
n_start               | int/null    | null    | Slice start index
keep_all              | bool        | false   | Disable all pruning
strict                | bool        | false   | Fail if all conformers removed
