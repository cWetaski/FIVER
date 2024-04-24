# FIVER
FIVER (FInite VolumE and Ray tracer) is a MATLAB-based set of tools which can be used to solve heat transfer problems involving radiation.

Radiation is resolved using a voxel-based ray tracer, which can (currently) handle opaque surfaces, nonscattering participating media, spectrally varying properties and refractive interfaces. 

In order to solve coupled heat transfer problems, the ray tracer has been coupled with Finite Volume Tool (FVTool) [1], the code for which is included in full in the subdirectory src/toolboxes/ (this is permitted under FVTool's BSD-2 license). Coupling with conduction has already been implemented for transient and equilibrium problems. But following the same methodology, it should not be too difficult to further couple the ray tracer with convection.

[1] Eftekhari, A.A. et al. (2015). FVTool: a finite volume toolbox for Matlab. Zenodo. http://doi.org/10.5281/zenodo.32745
