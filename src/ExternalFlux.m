%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07

classdef ExternalFlux < Flux
    %EXTERNALFLUX Defines an external flux which is applied to a voxel space
    %   The flux is defined positionally in 2 dimensions x and y (z-component of position from ray_gen function is ignored) 
    %   with direction vectors always having z-component >0, and then the flux is rotated according to the specified 
    %   boundary.
    %
    
    properties
        boundary; % (scalar String) ["left","right", "bottom","top","back","front"]: the boundary which the flux is applied to.
                  %         (left/right ==  -/+ x-axis == dim 1, bottom/top == -/+ y-axis == dim 2, back/front == -/+ z-axis == dim 3)
        end
    
    methods
        function [obj,voxel_space] = ExternalFlux(voxel_space,power,ray_gen_function,boundary)
            %EXTERNALFLUX Constructs an ExternalFlux Object and assigns it to the VoxelSpace
            obj.voxel_space = voxel_space;
            obj.power = power;
            obj.ray_gen_function = ray_gen_function;
            obj.boundary = boundary;
            voxel_space.fluxes = [voxel_space.fluxes,obj];
        end
        
        function [rays_pos,rays_dir] = GenerateRays(obj,N_rays)
            %   Detailed explanation goes here


            rays = obj.ray_gen_function(N_rays);

            if ~(all(rays(:,3)==0))
                disp('Warning: External flux ray gen function not defined with z-component of position == zero')
            end
            if any(rays(:,6)<0)
                disp('Warning: External flux ray gen function not defined with z-component of direction >= zero')
            end

            rays_pos = rays(:,1:2); % Take only x and y position components
            rays_dir = rays(:,4:6);

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

