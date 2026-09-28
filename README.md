# Benchmarking strain-level profiling of _Escherichia coli_ in short-read gut metagenomes (Galbraith _et al_., 2026)

## Summary

This repository contains the scripts from the Data Summary used for the benchmarking of short-read strain-level profiling tools (PanTax, PathoScope, StrainGE, Strainify, StrainR2 and StrainScan) on _Escherichia coli_. Reproducibility requirements for each tool are also provided. Supplementary data is available on the associated FigShare link from the paper (https://doi.org/10.6084/m9.figshare.32125474).

## Structure

- `envs/` - human-readable conda environment files used to define each tool’s dependencies. This includes InsilicoSeq, used for generating the simulated gut metagenomes and real-world spike-ins datasets.
- `locks/` - fully reproducible conda-lock files for exact environment recreation.
  - Example usage:
    - `conda-lock install -n <env_name> locks/<env_name>_conda_lock.yml`
    - `conda activate <env_name>`
  - Notes:
    -  PanTax requires Gurobi (v11). Gurobi is included in the lock file, but users must obtain and configure an appropriate Gurobi licence separately (see: https://www.gurobi.com/academics).
- `scripts/`
  - `InSilicoSeq_simulated_metagenomes/` - contains scripts and associated files (including abundance files) for simulating metagenomes using InSilicoSeq.
  - `R_visualizations_and_stats/` - contains scripts to recreate R visualizations and any associated statistical analysis used for figures in the manuscript. R scripts are grouped by results narrative - e.g., Fig. 1 and Fig. S1 are presented together. See FigShare repository for the required input data files.
  - `SLURM_jobs_simulated_metagenomes/` - contains per-tool scripts from the simulated metagenomes for job submissions on SLURM clusters.
    - The PanTax script is a file-driven batch processing loop due to Gurobi network requirements.
    - PathoScope has two scripts, one per core module used (MAP & ID), which should be run sequentially.
    - StrainGE has an additional extraction/normalisation script to process the results.
  - `real_world_spike_ins_metagenomes/` - contains per-tool scripts from the real-world spike-ins final benchmarking dataset for DB building, metagenome profiling and computational benchmarking.

## Preprint

Now available on BioRxiv, see https://doi.org/10.64898/2026.05.19.726160
