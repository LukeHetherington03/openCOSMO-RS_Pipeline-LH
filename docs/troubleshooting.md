# Troubleshooting Guide

openCOSMO-RS Pipeline

---

## How to diagnose a failed run

**Step 1 — Check the pipeline state:**
```bash
pl r status <request-id>
```
This shows which stage failed and the top-level error message.

**Step 2 — Read the stage log:**
```bash
pl r logs <request-id>
```
Look for `[FAILED]`, `[ERROR]`, and `[STAGE] ... FAILED` lines. The `[FAILED]` lines include an `→ error:` continuation row with the specific error.

**Step 3 — Check the pipeline warnings summary:**
```
pipeline_data/requests/<request-id>/final_outputs/pipeline_warnings.txt
```
Or in the request log:
```
pipeline_data/requests/<request-id>/request.log
```

**Step 4 — Examine the raw job state:**
```
pipeline_data/requests/<request-id>/jobs/<job-id>/job_state.json
```
Shows item-level counts (completed, failed, pending) and the stored error message.

---

## Path truncation errors (XTB / gXTB)

**Symptoms:**
- XTB or gXTB fails with a file not found error
- Error message contains truncated-looking paths (e.g. path ends mid-directory-name)
- Works fine when software is in a shorter path

**Cause:**
XTB and gXTB are written in Fortran and use fixed-length character buffers for file paths. If any path (to the executable, the working directory, or the input files) exceeds the buffer length, it is silently truncated and the tool fails.

**Fix:**
Move the software to a shorter path. The recommended location is `~/software/`:
```
~/software/xtb/bin/xtb
~/software/g-xtb-main/binary/gxtb
~/software/orca_6_1_1/orca
~/software/crest/crest
```
Update `config/paths.json` with the new paths. Do **not** install software inside the project directory tree — the combined depth can easily exceed the buffer limit.

---

## ORCA segfaults or silent failures

**Symptoms:**
- ORCA jobs fail immediately with exit code 139 (segfault) or no output
- `[FAILED]` log shows `returncode != 0` with no useful error message
- Works with shorter paths or fewer cores

**Possible causes and fixes:**

| Cause | Fix |
|-------|-----|
| Path too long (similar to XTB issue) | Move ORCA to `~/software/orca_6_1_1/` |
| `maxcore` too high — ORCA exceeds available RAM | Reduce `maxcore` in `orcacosmo_defaults.json` (default: 2000 MB) |
| `cores_per_item` × `maxcore` exceeds node RAM | Reduce `cores_per_item` or `maxcore` |
| Scratch directory on a slow/full filesystem | Check disk space; ensure `base_dir` points to a fast disk |
| Wrong ORCA version | Verify the executable works: `~/software/orca_6_1_1/orca --version` |

---

## Missing COSMO files / solubility stage fails

**Symptoms:**
- Solubility stage fails with "no conformers" or "COSMO file not found"
- `orcacosmo_summary.json` exists but some entries have `status: "failed"`

**Checks:**
1. Confirm the orcacosmo stage completed successfully: `pl r status <id>` should show `orcacosmo` as completed.
2. Check `final_outputs/pipeline_warnings.txt` for orcacosmo warnings — fallback triggers indicate basis set issues.
3. Verify the COSMO output files exist:
   ```
   pipeline_data/requests/<id>/jobs/<orcacosmo-job-id>/outputs/orcacosmo_outputs/
   ```
4. To add a fallback basis, add a second `orcacosmo` stage in your pipeline spec with the desired basis set and `"fallback_only": true`. That stage will only run ORCA on conformers that do not yet have a valid `.orcacosmo` file.

---

## Solvent not found

**Symptoms:**
- Solubility stage skips one or more combos with a warning
- `predicted_solubility: null` for some solvents

**Cause:**
The solvent name in the solvent list JSON does not match a directory under `CONSTANT_FILES/solvents/`.

**Fix:**
Check the available solvents:
```bash
ls CONSTANT_FILES/solvents/
```
The name in the solvent list must exactly match the directory name (case-sensitive). See [solvent_list_format.md](solvent_list_format.md) for how to add new solvents.

---

## Worker not processing requests

**Symptoms:**
- `pl r list` shows requests as `pending` but they never run
- `pl q status` shows the worker is not running

**Fix:**
```bash
pl q start
```

If the worker appears to be running but requests are stuck:
```bash
pl q status    # check for errors
pl q stop
pl q start
```

If the worker process died unexpectedly, its lock file may be stale. Check and clear:
```bash
pl q status    # will report if worker is unresponsive
```

---

## Request fails immediately after submission

**Symptoms:**
- Stage fails before processing any items
- Error refers to a missing file or import error

**Common causes:**

| Error | Cause | Fix |
|-------|-------|-----|
| `ModuleNotFoundError` | Python dependency missing | `pip install -r requirements.txt` |
| `stage_input not found` | Previous stage output missing | Check the previous stage completed; look at job outputs |
| `config/paths.json not found` | Not running from project root | `cd` to the project root before running `pl` |
| `executable not found` | Path in `paths.json` is wrong | Run `pl env software` to validate all executables |

---

## Code changes not taking effect

**Symptoms:**
- You edited a source file but the pipeline still runs the old code
- Bug appears fixed in the file but still occurs at runtime

**Cause:**
The queue worker is a long-running process. Python loads modules once and caches them in `__pycache__/`. Changes to `.py` files are not picked up until the worker restarts.

**Fix:**
```bash
pl q stop    # wait for running items to finish checkpoints
pl q start
```

Or for an immediate restart:
```bash
pl q kill
pl q start
```

---

## Checkpoint corruption

**Symptoms:**
- Stage fails with a JSON parse error on a checkpoint file
- Resume skips items that should not be skipped

**Cause:**
A checkpoint was partially written during a crash (interrupted between open and write). The pipeline uses atomic writes to minimise this risk, but a full system crash can still cause corruption.

**Fix:**
Delete the corrupted checkpoint file:
```
pipeline_data/requests/<id>/jobs/<job-id>/outputs/checkpoints/<item_id>.json
```
Then resubmit the request — the item will be re-run from scratch.

To delete all checkpoints for a job and re-run the entire stage, delete the `checkpoints/` directory.

---

## Environment validation

Run the full environment check at any time:
```bash
pl env check
```

This validates:
- All executables configured in `paths.json` exist and are executable
- openCOSMO-RS Python source and C++ bindings are accessible
- CONSTANT_FILES directories exist
- Resource allocation configuration is valid

For targeted checks:
```bash
pl env software      # executables only
pl env resources     # openCOSMO paths and constant files
```

---

## Useful paths for manual inspection

```
config/paths.json                               global software paths
config/resource_allocation.json                 core budget configuration
pipeline_data/requests/                         all request data
pipeline_data/requests/<id>/request.log         high-level pipeline log
pipeline_data/requests/<id>/request_state.json  pipeline state
pipeline_data/requests/<id>/jobs/<job-id>/
    job_state.json                              job status and item counts
    stage_logs/stage.log                        full stage log
    stage_logs/stage_context.log                configuration snapshot at startup
    inputs/parameters.json                      all resolved parameters for this job
    outputs/checkpoints/                        per-item checkpoints
    outputs/                                    stage outputs
pipeline_data/requests/<id>/final_outputs/      results and warnings summary
```
