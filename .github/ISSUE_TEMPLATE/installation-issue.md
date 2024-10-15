---
name: Installation Issue
about: Trouble installing WarpX? Let us help you.
title: ''
labels: [install, question]
assignees: ''

---

**Where are you trying to install and run WarpX?**
A clear and concise description of where you like to run.

**Software Environment**
- version of WarpX: [e.g., latest or 24.09]
- operating system: [name and version]
- machine: [Are you running on a public cluster? It's likely we compute on it as well!]

**What have you tried so far?**
- installing WarpX via: [conda-forge, spack, pip, brew, from source, module system, ...]

**Installation from Source**
Are you installing form source? If so, please provide:
- the output of `cmake --fresh -S . -B build`
- the output of `cmake --build build`

Let us also know more about your software environment in that case:
- version of CMake: [e.g. 3.24.0]
- name and version of your C++ compiler: [e.g., GNU 11.3 with NVCC 12.0.76]
- name and version of Python implementation: [e.g. CPython 3.12]
- version of ADIOS2: [e.g. 2.10.0]
- version of FFTW: [e.g. 3.3.10]
- version of HDF5: [e.g. 1.14.0]
- name and version of MPI: [e.g. OpenMPI 4.1.1]
- name of other relevant dependency names and version: [e.g. NumPy 2.0, Ascent 0.8.0, etc.]

**Screenshots**
If applicable, add screenshots to help explain your problem.
