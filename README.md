# Benchmarking strain-level profiling of _Escherichia coli_ in short-read gut metagenomes

## Summary

This repository contains the scripts from the Data Summary used for benchmarking short-read strain-level profiling tools.

## Structure

- `scripts/InSilicoSeq_simulated_metagenomes` - contains scripts and associated files for simulating metagenomes using InSilicoSeq. The simulated metagenomes generated from these scripts are also provided on FigShare.
- `scripts/SLURM_jobs_simulated_metagenomes` - contains per-tool scripts from the simulated metagenomes for job submissions on SLURM clusters. Note PanTax has a script included here for completeness but this is not a SLURM job script, rather a file-driven batch processing loop. Note that PathoScope has two scripts, one per module (MAP & ID), which should be run sequentially.
