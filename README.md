# MHIT36

Code for direct numerical simulation of Navier-Stokes equation coupled with a phase-field method (ACDI) for interface description.

~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████       
████  ████ ██   ██ ██    ██         ██ ██              
██ ████ ██ ███████ ██    ██     █████  ███████   
██  ██  ██ ██   ██ ██    ██         ██ ██    ██     
██      ██ ██   ██ ██    ██    ██████   ██████        
~~~

Multi-GPU version of MHIT36 using cuDecomp (Nvidia only)

If you use this code, please cite the following work: 
```bibtex
  @article{roccon2025,
  title   = {MHIT36: A Phase-Field Code for Gpu Simulations of Multiphase Homogeneous Isotropic Turbulence},
  author  = {Roccon, A. and Enzenberger, L. and Zaza, D. and Soldati, A.},
  journal = {SSRN},
  year    = {2025},
  doi     = {http://dx.doi.org/10.2139/ssrn.5264052}
}
```

![Test](val/intro.png)


## Check list of features implemented in MHIT36

- Poisson solver (transposition + halo update) ✅
- Poisson solver validation (periodic solutions) ✅
- Read input files ✅
- Skeleton of the code  ✅
- Halo updates test with CUDA ✅
- Poisson solver scaling ✅
- Halo updates test with host_data use_device ✅
- Flow field initialization ✅
- Phase-field initialization ✅
- Phase-field method (ACDI) ✅
- Forcing ✅
- HIT validation ✅
- Full code scaling ✅
- MPI writing (no halo)  ✅
- MPI reading (no halo)  ✅
- Umax via MPI reduction ✅
- Surface tension forces ✅
- Remove mean flow via MPI all reduce ✅

## Run the code

- Compile first the cuDecomp library using *_lib.sh, the resulting modules and library will be located in cuDecomp/build/lib and cuDecomp/build/include
- Double check cuDecomp building is fine (must be compiled using HPC-SDK)
- Folder multi: contains the source-code of the multi GPU version of the code. Use local.sh, leo.sh or mn5.sh to compile and run the code; the multi GPU version relies on cuDecomp for pencils transpositions and halo exchanges.
- Autotuning of the multi-GPU version: Default pr=0 and pc=0 enables autotuging (when cuDecomp is initialized), cuDecomp will perform an autotuning at the start finding the best decomposition (the only input is the total number of tasks). In this way, everyhting is automatic and the code does not need to be recompiled when changing the number of MPI processes.
- A conditional compilation flag is used to enable or not the phase-field module. By default is single-phase only.


## Reference performance

Performance (NS only)
* 256 x 256 x 256 | 4 x A100@Leonardo  |   16 ms/timestep
* 512 x 512 x 512 | 4 x A100@Leonardo  |  270 ms/timestep 
* 1024 x 1024 x 1024 | 32 x A100@Leonardo | 357 ms/timestep 
* 2048 x 2048 x 2048 | 128 x A100@Leonardo | 820 ms/timestep
* 256 x 256 x 256 | 4 x H100@MN5-ACC   |   13 ms/timestep
* 512 x 512 x 512 | 4 x H100@MN5-ACC   |  230 ms/timestep
* 1024 x 1024 x 1024 | 32 x H100@MN5-ACC  | 320 ms/timestep

* 2048 x 2048 x 2048 | 512 x H100@MN5-ACC  | 259 ms/timestep

Phase-field introduces about 15% of overhead compared to NS only.

## Scaling

Strong scaling results obtained on Leonardo (4 x A100 64 GB x node) and MN5 (4 x H100 40 GB x node)
* Tested from 1 node up to 128 nodes (Leonardo)
* Tested from 1 node up to 256 nodes (MN5-ACC)
* Grid from 64 x 64 x 64 up to 2048 x 2048 x 2048
* Very similar scaling for both NS and NS+ACDI

![Scal](val/scaling.png)


## Validation

Benchamrk present in "W.M.VanRees, A.Leonard, D.Pullin, P.Koumoutsakos, A comparison of vortex and pseudo-spectral methods for the simulation of periodic vortical flows at high Reynolds numbers,J. Comput. Phys.2 30(8)(2011)2794–2805" and also Used in CaNS.

Time evolution of the viscous dissipation:

![Test](val/val.png)


## Contributing

We welcome all contributions that can enhance MHIT36, including bug fixes, performance improvements, and new features. 
If you would like to contribute, please contact aroccon or open an Issue in the repository.