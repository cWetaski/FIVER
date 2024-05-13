classdef ExternalFlux
    %EXTERNALFLUX Defines an external flux which is applied to a voxel space
    %   The flux is defined in 2 Dimensions x and y, with direction vectors always having z-component >0,
    %   and then the flux is rotated according to the actual boundary which is applied to.
    %   It is assumed that the arguments to the Inverse PDF functions are uniformly generated random numbers in [0,1]
    %   It is also assumed that the Inverse PDF functions are vectorized (for generating many rays)
    
    properties
        voxel_space; % VoxelSpace object which the external flux is applied to
        boundary; % (scalar String) ["left","right", "bottom","top","back","front"]: the boundary which the flux is applied to.
                  %         (left/right ==  -/+ x-axis == dim 1, bottom/top == -/+ y-axis == dim 2, back/front == -/+ z-axis == dim 3)
        power; % (scalar) [W]: The total power of the external flux.  
        ray_gen_function; % (Function of scalar N_rays): Inverse cumulative distribution function for generating a ray originating in the x,y plane
                             %                              The generated rays have components [x,y, dx,dy,dz] with dz > 0
        end
    
    methods
        function [obj,voxel_space] = ExternalFlux(voxel_space,boundary,power,ray_gen_function)
            %EXTERNALFLUX Constructs an ExternalFlux Object and assigns it to the VoxelSpace
            obj.voxel_space = voxel_space;
            obj.boundary = boundary;
            obj.power = power;
            obj.ray_gen_function = ray_gen_function;
            voxel_space.external_fluxes = [voxel_space.external_fluxes,obj];
        end
        
        function [rays_pos,rays_dir] = GenerateRays(obj,N_rays)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            rays_xy = obj.ray_gen_function(N_rays);
            rays_pos = rays_xy(:,1:2);
            rays_dir = rays_xy(:,3:5);

            %% Rotate the generated points and vectors around the center point of the domain so that they are in the appropriate plane
            % https://math.stackexchange.com/questions/2093314/rotation-matrix-of-rotation-around-a-point-other-than-the-origin
            size_VS = obj.voxel_space.size;
            center_point = size_VS/2; % Minimum of domain is always 0 in all 3 axes
            switch obj.boundary % Define rotation matrix based on boundary that flux is entering through
                case 'left' % Rotate so that points lie in y-z plane, positive z-component of direction becomes positive x-component
                            % i.e. counterclockwise rotation of 90 degrees about the y-axis
                    M_rot = [ 0,  0,  1;...
                              0,  1,  0;...
                             -1,  0,  0]; 

                case 'right' % Rotate so that points lie in y-z plane, positive z-component of direction becomes negative x-component
                             % i.e. clockwise rotation of 90 degrees about the y-axis
                    M_rot = [ 0,  0, -1;...
                              0,  1,  0;...
                              1,  0,  0];

                case 'bottom' % Rotate so that points lie in x-z plane, positive z-component of direction becomes positive y-component 
                              % i.e. clockwise rotation of 90 degrees about the x-axis
                    M_rot = [ 1,  0,  0;...
                              0,  0,  1;...
                              0, -1,  0];

                case 'top' % Rotate so that points lie in x-z plane, positive z-component of direction becomes negative y-component 
                           % i.e. counterclockwise rotation of 90 degrees about the x-axis
                    M_rot = [ 1,  0,  0;...
                              0,  0, -1;...
                              0,  1,  0];

                case 'back' % No rotation required, this is the x-y plane with positive z-component pointing into the domain
                    M_rot = eye(3); % Identity matrix

                case 'front' % Rotate 180 degrees around either the x or y axis (or both)
                             % we will rotate about the y axis since it's easier to rotate my right-hand in this way when I'm doing the right-hand rule
                    M_rot = [-1,  0,  0;...
                              0,  1,  0;...
                              0,  0, -1];
            end
            rays_pos = [rays_pos,zeros(N_rays,1)];
            rays_pos = (M_rot*(rays_pos - center_point)')'+center_point; 
            rays_dir = (M_rot*rays_dir')';
        end
    end
end

