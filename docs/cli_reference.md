# Command‑Line Interface (CLI) Reference  
openCOSMO‑RS Pipeline

This document provides a complete reference for the `pl` command‑line interface.  
The CLI is used for:

- Managing the execution queue  
- Inspecting Requests  
- Validating the environment  

The CLI **does not** create or run Requests.  
Requests are created via the Python entrypoints:

```
python3 main.py
python3 main_resume.py
python3 main_continue.py
```

---

# 1. Setup

Add the alias to your `.bashrc`:

```bash
alias pl='python3 -m modules.cli.cli'
```

Reload:

```bash
source ~/.bashrc
```

You can now run:

```
pl <command>
```

---

# 2. Command Groups

The CLI has three command groups:

| Group | Prefix | Purpose |
|-------|--------|----------|
| Queue commands | `pl q` | Manage the execution queue & worker |
| Request commands | `pl r` | Inspect and annotate Requests |
| Environment commands | `pl env` | Validate installation & dependencies |

Each group is documented below.

---

# 3. Queue Commands (`pl q`)

Queue commands manage queued Requests and the worker process.

### **Start the worker**
```
pl q start
```
Starts the queue worker.  
The worker continuously processes queued Requests.

### **Stop the worker**
```
pl q stop
```
Stops the worker gracefully.

### **Check queue status**
```
pl q status
```
Shows:
- Worker status  
- Pending Requests  
- Running Requests  
- Completed Requests  

### **List queued Requests**
```
pl q list
```
Displays all Requests currently in the queue.

### **Cancel a queued Request**
```
pl q cancel <id>
```
Removes a Request from the queue.

### **Reprioritise a queued Request**
```
pl q reprio <id> <prio>
```
Changes the priority of a queued Request.

---

# 4. Request Commands (`pl r`)

Request commands inspect and annotate Requests.  
They do **not** execute or modify scientific results.

### **List all Requests**
```
pl r list
```

### **Check Request status**
```
pl r status <id>
```

### **Show Request metadata**
```
pl r info <id>
```

### **View Request logs**
```
pl r logs <id>
```

### **Add a note to a Request**
```
pl r note <id> "text"
```

### **Tag a Request**
```
pl r tag <id> tag1 tag2
```

### **Pin a Request**
```
pl r pin <id>
```
Marks a Request as important.

### **Publish a Request**
```
pl r publish <id>
```
Marks a Request as ready for sharing.

### **Archive a Request**
```
pl r archive <id>
```
Moves a Request to the archive.

### **Trash a Request**
```
pl r trash <id>
```
Marks a Request for deletion.

---

# 5. Environment Commands (`pl env`)

Environment commands validate your installation.

### **Full environment check**
```
pl env check
```

### **Validate external executables**
```
pl env software
```
Checks ORCA, XTB, gXTB, CREST, etc.

### **Validate openCOSMO paths + constant files**
```
pl env resources
```

### **Validate chemistry JSON files**
```
pl env chemistry
```

### **Validate required Python packages**
```
pl env pip
```
Example output:

```
======================================================================
Pip Package Validation
======================================================================
✓ numpy                     version 1.26.4
✓ scipy                     version 1.14.1
✓ pandas                    version 2.3.3
✓ rdkit                     version 2025.09.1
✓ openbabel                 version 3.1.0
```

### **Pretty table output**
```
pl env table
```

---

# 6. Summary

The CLI provides:

- Queue management (`pl q …`)
- Request inspection (`pl r …`)
- Environment validation (`pl env …`)

Requests are always created via:

```
python3 main.py
python3 main_resume.py
python3 main_continue.py
```

The CLI never runs or creates Requests — it manages them.

---

# 7. Related Documents

- `docs/execution_and_queueing.md`  
- `docs/installation.md`  
- `docs/configuration.md`  
- `docs/pipeline_architecture.md`  

