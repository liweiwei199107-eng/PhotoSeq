# Windows R and Python environments

The downstream analysis environment was:

- Windows 11 x64, build 26200
- R 4.6.0 (2026-04-24, UCRT)
- Bioconductor 3.23
- UTF-8 system code page

The project-level `../renv.lock` records 300 installed dependency packages. The active `msigdbr` version is 26.1.1. `PhotoSeq_R_package_versions.csv` is the corresponding non-duplicated summary of the principal analysis packages, and `R_sessionInfo.txt` records the R platform.

`required_R_packages.csv` lists every non-base package called directly by the analysis scripts. `renv.lock` provides the complete dependency record.

## Restore the recorded environment

Install R 4.6.0, open a terminal in `02_downstream_analysis/` and run:

```powershell
R -e "install.packages('renv', repos='https://cloud.r-project.org')"
R -e "renv::restore()"
```

Restoring the lock file requires HTTPS access to the recorded CRAN and Bioconductor-compatible repositories.

The equivalent commands from an R session are:

```r
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore()
```

After restoration, verify the packages called directly by the scripts:

```powershell
R -e "renv::load(); source('environment/check_required_packages.R')"
```

`run_all_analysis.R` activates the restored project library and passes it to each numbered child process automatically.

## Python environment

The recorded Python environment used Python 3.12.14 and pandas 3.0.1. From the
`02_downstream_analysis/` directory, install the recorded dependency with:

```powershell
python -m pip install -r environment\python_requirements.txt
```

The PPI network-scoring helper itself uses only the Python standard library. A
live STRING query requires HTTPS access to
`https://string-db.org`; the default reproducibility run uses the frozen
response under `reference/ppi/` and does not require network access.
