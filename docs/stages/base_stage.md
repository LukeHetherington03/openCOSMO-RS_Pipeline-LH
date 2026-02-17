# BaseStage  
openCOSMO‑RS Pipeline

All scientific stages in the pipeline inherit from `BaseStage`.  
This class provides a unified lifecycle, logging, item tracking, strict‑mode helpers, and canonical input/output handling.

Understanding `BaseStage` is essential for developing or modifying any stage.

---

# 1. Purpose

`BaseStage` provides:

- A consistent execution lifecycle  
- Automatic logging  
- Automatic job state updates  
- Canonical input/output helpers  
- Strict‑mode enforcement  
- Item‑level progress tracking  
- Error handling and pipeline state updates  

Every stage in `modules/stages/` subclasses `BaseStage` and implements:

```python
def execute(self):
    ...
```

---

# 2. Lifecycle

Stages are executed through:

```python
stage.run()
```

This wrapper:

1. Logs stage start  
2. Marks the job as running  
3. Calls `execute()`  
4. Marks the job as complete  
5. Updates pipeline_state.json  
6. Handles exceptions  
7. Marks job as failed if needed  

This ensures consistent behaviour across all stages.

---

# 3. Canonical Input & Output

### Input
Stages may declare a required input:

```python
stage_input = self.parameters.get("stage_input")
```

If missing:

```python
self.fail("Stage requires 'stage_input' but none was provided.")
```

### Output
Stages declare their canonical output using:

```python
self.set_stage_output("cleaned.csv")
```

This sets:

```
outputs/<filename>
```

and registers the canonical output for the pipeline.

---

# 4. Logging

`BaseStage` provides:

```python
self.log(msg)
self.log_header(title)
self.log_section(title)
```

All logs are written to:

```
jobs/J-.../stage.log
```

Logging is consistent across all stages.

---

# 5. Item Tracking

Stages that operate on multiple items (e.g., molecules) use:

```python
self.set_items(list_of_items)
self.update_progress(item)
```

This updates:

```
job_state.json
```

with:

- pending_items  
- completed_items  
- failed_items  

This enables deterministic resume.

---

# 6. Strict Mode

Stages may enforce strict validation:

```python
if self.strict("cleaning"):
    ...
```

Strict mode is configured in:

```
config/paths.json
```

Example:

```json
"cleaning": { "strict": true }
```

---

# 7. Path Helpers

Stages use:

```python
self.input_path("file")
self.output_path("file")
```

These resolve to:

```
jobs/J-.../inputs/file
jobs/J-.../outputs/file
```

This ensures deterministic file layout.

---

# 8. Error Handling

To fail a stage:

```python
self.fail("Something went wrong")
```

This:

- Logs the error  
- Marks the job as failed  
- Updates pipeline_state.json  
- Raises RuntimeError  

---

# 9. Developer Template

A minimal stage looks like:

```python
from modules.stages.base_stage import BaseStage

class ExampleStage(BaseStage):

    def execute(self):
        self.log_section("Loading input")
        input_file = self.get_stage_input()

        self.log_section("Processing")
        # ... do work ...

        self.log_section("Writing output")
        output_path = self.set_stage_output("example_output.json")
        with open(output_path, "w") as f:
            json.dump({"ok": True}, f)
```

---

# 10. Summary

`BaseStage` ensures:

- Consistent lifecycle  
- Deterministic behaviour  
- Reproducibility  
- Safe resume  
- Clean logging  
- Canonical outputs  

Every stage in the pipeline builds on this foundation.

