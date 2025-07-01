# LSF_utils

This directory contains helper scripts and job control wrappers used to manage LSF-based job scheduling in the somatic mutation detection pipeline.

These utilities are primarily used by scripts in `Mu2_1stCall/` and `Mu2_2ndCall/` to automate parallel job submission, file lookup by sample tags, and job dependency tracking.

---

## 🔁 Job Submission Wrappers

- **`BsubS`**  
  Standard job submission wrapper for quick tasks (e.g., indexing)

- **`BsubM`**  
  Medium or memory-intensive job submission wrapper

- **`BsubL`**  
  Long-duration or high-memory job submission wrapper

> These wrappers encapsulate `bsub` with preconfigured resource parameters (e.g., `-n`, `-M`, and queue).

---

## 🔍 File Tag Utilities

- **`File_pickByTag.sh`**  
  Retrieves file path(s) from a list using a sample tag. Returns matched path(s) if exactly one hit is found.

- **`File_pickByTag_onlyExitStatus.sh`**  
  Same as above, but returns only exit code (0 for success, 1 for failure). Used for logical gating in shell scripts.

---

## ⏳ Job Dependency Tools

- **`WaitSignalMaker.sh`**  
  Converts a job ID list (e.g., from `bsub`) into a properly formatted LSF `-w` dependency string (e.g., `done(job1) && done(job2)`).

---

## Example Usage

```bash
# Retrieve BAM path for a sample tag
File_pickByTag.sh BamList.txt Sample001

# Check if tag uniquely matches one path
File_pickByTag_onlyExitStatus.sh BamList.txt Sample001

# Create wait signal for downstream dependent job
WaitSignalMaker.sh JobID_Mu2Call_PID01.txt
```

---

## Requirements

- LSF Job Scheduler (`bsub`)
- Bash 4+
- File lists formatted with one record per line
