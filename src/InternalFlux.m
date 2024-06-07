classdef InternalFlux
    %INTERNALFLUX Defines an internal flux within a voxel space
    %   Defines a position and direction distribution for a source flux within a voxel space as well as the overall
    %   power (in W) associated with the flux
    
    properties
        voxel_space;      % VoxelSpace object which the external flux is applied to
        power;            % (scalar) [W]: The total power of the external flux.  
        ray_gen_function; % (Function of scalar N_rays): Generates N_rays according to some probability distribution
                          %                              The generated rays have components [x,y,z, dx,dy,dz]
        end
    
    methods
        function [obj,voxel_space] = InternalFlux(voxel_space,power,ray_gen_function)
            %EXTERNALFLUX Constructs an ExternalFlux Object and assigns it to the VoxelSpace
            obj.voxel_space = voxel_space;
            obj.power = power;
            obj.ray_gen_function = ray_gen_function;
            voxel_space.fluxes = [voxel_space.fluxes,obj];
        end
        
        function [rays_pos,rays_dir] = GenerateRays(obj,N_rays)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            rays = obj.ray_gen_function(N_rays);
            rays_pos = rays(:,1:3);
            rays_dir = rays(:,4:6);
        end
    end
end

