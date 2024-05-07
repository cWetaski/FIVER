# FIVER
FIVER (FInite VolumE and Ray tracer) is a MATLAB-based set of tools which can be used to solve heat transfer problems involving radiation. This code is an extension of the work of Sebastian Sas-Brunser, described in [1].

Radiation is resolved using a voxel-based ray tracer, which can (currently) handle opaque surfaces, nonscattering participating media, spectrally varying properties and refractive interfaces. 

In order to solve coupled heat transfer problems, the ray tracer has been coupled with Finite Volume Tool (FVTool) [2], the code of which is included in full in the subdirectory src/toolboxes/. Coupling with conduction has already been implemented for transient and equilibrium problems. But following the same methodology, it should not be straightforward to further couple the ray tracer with convection. 

[1] Sas-Brunser, S. and Steinfeld, A. (2023). Design and Optimization of Hierarchically Ordered Porous Structures for Solar Thermochemical Fuel Production Using a Voxel-Based Monte Carlo Ray-Tracing Algorithm, ACS Engineering Au 2023 3 (5), 326-334, DOI: 10.1021/acsengineeringau.3c00013

[2] Eftekhari, A.A. et al. (2015). FVTool: a finite volume toolbox for Matlab. Zenodo. http://doi.org/10.5281/zenodo.32745

## How to use the code
To add all the functions to your MATLAB path, simply open the FIVER.proj file (Alternatively you can add all the files in the /src/ folder to your path manually, and install FVTool as a MATLAB add-on).

Several example codes are stored in /test/, each function also has detailed descriptions. 

TODO: Add more example codes

TODO: Wiki documenting overall code methodology. 
