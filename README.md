# TCF36

Extended version of MHIT36 for turbulent channel flow.
Finite-difference code based on the fractional step method.
Solution of Navier-Stokes equations + phase-field method (ACDI).

Time integration, two options:
- All the terms explicit using RK3 (for both NS + PFM, working)
- CN for z-diffusive terms and RK3 for the rest of the terms (to be implemented, only skeleton is present), PFM is explicit.


If you use this code, please cite the following work: 
```bibtex
  @article{roccon2025,
  title   = {MHIT36: A Phase-Field Code for Gpu Simulations of Multiphase Homogeneous Isotropic Turbulence},
  author  = {Roccon, A. and Enzenberger, L. and Zaza, D. and Soldati, A.},
  journal = {Computer Physics Communications (in press)},
  year    = {2025},
  doi     = {https://doi.org/10.1016/j.cpc.2025.109804}
}
```


## Check list of TCF36
- Boundary condition for no-slip at the two walls ✅
- Laminar solution (no need of TDMA) ✅
- TDMA ✅
- TDMA validation ✅ 
- Turbulent channel flow ✅ 
- Stretched grids 🚧
- Implicit diffusion along z (skeleton and flag introduced) 🚧

# How to run the code
## 1. Compile cuDecomp
- First, build the **cuDecomp** library using the corresponding `*_lib.sh` script.  
- The resulting files will be located in:
  - **Library:** `cuDecomp/build/lib`  
  - **Headers:** `cuDecomp/build/include`  
- ⚠️ Ensure that cuDecomp is compiled using **NVIDIA HPC-SDK**.  

## 2. Multi-GPU Source Code
- The folder **`multi/`** contains the source code for the **multi-GPU** version.  
- To compile and run:
  - Use `local.sh`, `leo.sh`, or `mn5.sh` depending on your system.  
- The multi-GPU version relies on **cuDecomp** for:
  - Pencil decompositions and transpositions  
  - Halo exchanges between GPUs  

## 3. Autotuning of Multi-GPU Decomposition
- By default, `pr=0` and `pc=0` enable **autotuning** when cuDecomp is initialized.  
- cuDecomp will automatically determine the best process decomposition at runtime.  
- The only required input is the **total number of MPI tasks**.  
- ✅ This means you **do not need to recompile** the code when changing the number of MPI processes.  

## 4. Conditional Compilation Flags
- **Phase-field module**: Can be enabled or disabled. By default, only single-phase is used.  
- **Implicit diffusion along z**: Can be enabled or disabled. ⚠️ This feature is **not yet implemented**.


## Turbulent channel flow 

![Test](val/tcf.png)

## Performance and resolution tested (NS only)

- 256 x 128 x 200 - 31 ms/iter - 2 x RTX5000 16GB 
- 2048 x 768 x 576 - 323 ms/iter - 4 x A100 64 GB 

## Contributing

We welcome all contributions that can enhance TCF36, including bug fixes, performance improvements, and new features. 
If you would like to contribute, please contact me or open an Issue in the repository.