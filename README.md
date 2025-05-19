# EvoRSbent

Code for the poster "Evolving Cryptographic Boolean Functions with Reaction Systems", accepted at Gecco 2025.

## Repository Structure

- `Bool_module.py`: contains functions related to Boolean operations and Boolean functions.
- `fitness.py`: containts different definition of fitness functions.
- `RSmodule.py`: implements the class of Reaction Systems.
- `EvoRSmodule_newmuta_nobal.py`: core module implementing the evolutionary algorithm.
- `EvoRS_run_fitnl_new_mutation_nobal.py`: script to execute the evolutionary algorithm with specified parameters.
- `run_tests_k=6_bent.sh`, `run_tests_k=8_bent.sh`, `run_tests_k=10_bent1.sh`: shell scripts to run tests for different values of `k`.
