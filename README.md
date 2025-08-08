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
- Turbulent channel flow 🚧 (looks promising)
- Stretched grids 🚧
- Implicit diffusion along z (skeleton and flag introduced) 🚧

## Run the code

- Compile first the cuDecomp library using *_lib.sh, the resulting modules and library will be located in cuDecomp/build/lib and cuDecomp/build/include
- Double check cuDecomp building is fine (must be compiled using HPC-SDK)
- Folder multi: contains the source-code of the multi GPU version of the code. Use local.sh, leo.sh or mn5.sh to compile and run the code; the multi GPU version relies on cuDecomp for pencils transpositions and halo exchanges.
- Autotuning of the multi-GPU version: Default pr=0 and pc=0 enables autotuging (when cuDecomp is initialized), cuDecomp will perform an autotuning at the start finding the best decomposition (the only input is the total number of tasks). In this way, everyhting is automatic and the code does not need to be recompiled when changing the number of MPI processes.
- A conditional compilation flag is used to enable or not the phase-field module. By default is single-phase only.
- A conditional compilation flag is used to enable or not implicit diffusion integration along z; this feature is not yet implemented.

## Turbulent channel flow 

![Test](val/tcf.png)

## Performance and resolution tested (NS only)

- 256 x 128 x 200 - 45 ms/iter - 2 x RTX5000 16GB 
- 2048 x 768 x 576 - 410 ms/iter - 4 x A100 64 GB 

## Contributing

We welcome all contributions that can enhance TCF36, including bug fixes, performance improvements, and new features. 
If you would like to contribute, please contact me or open an Issue in the repository.