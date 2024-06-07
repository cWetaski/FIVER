classdef (Abstract) Flux
    %FLUX Parent class of internal and external fluxes which can be assigned to a VoxelSpace object. 
    %   Fluxes are defined using:
    %       1) A ray generator function which generate ray positions and directions according to some underlying 
    %          probability distribution function
    %       2) An associated power (in W)
    
    properties
        voxel_space;      % VoxelSpace object which the external flux is applied to
        power;            % (scalar) [W]: The total power of the external flux.  
        ray_gen_function; % (Function of scalar N_rays): Generates N_rays according to some probability distribution
                          %                              The generated rays have components [x,y,z, dx,dy,dz]
                          %                              In the case of an ExternalFlux, the ray_gen_function should always give 0 z-component
    end
    
    methods (Abstract)
        GenerateRays(obj,N_rays);
    end
end

