% ©2024 ETH Zurich, Charles Wetaski, Sebastian Sas Brunser, Emiliano Casati
%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-12-10

classdef ExternalFlux < Flux
    %EXTERNALFLUX Defines an external flux which is applied to a voxel space
    %   The flux is defined positionally in 2 dimensions x and y (z-component of position from ray_gen function is ignored) 
    %   with direction vectors always having z-component >0, and then the flux is rotated according to the specified 
    %   boundary.
    %
    
    properties
        boundary; % (scalar String) ["left","right", "bottom","top","back","front"]: the boundary which the flux is applied to.
                  %         (left/right ==  -/+ x-axis == dim 1, bottom/top == -/+ y-axis == dim 2, back/front == -/+ z-axis == dim 3)
        size_VS;     % (1x3 double) size of the voxel space
    end
    
    methods
        function obj = ExternalFlux(power,ray_gen_function,boundary,size_VS)
            %EXTERNALFLUX Constructs an ExternalFlux Object
            obj.power = power;
            obj.ray_gen_function = ray_gen_function;
            obj.boundary = boundary;
            obj.size_VS = size_VS;
        end
        
        function [rays_pos,rays_dir] = GenerateRays(obj,N_rays)
            %   Generate rays
            rays = obj.ray_gen_function(N_rays);

            if ~(all(rays(:,3)==0))
                disp('Warning: External flux ray gen function not defined with z-component of position == zero')
            end
            if any(rays(:,6)<0)
                disp(['Warning: External flu' ...
                    'x ray gen function not defined with z-component of direction >= zero'])
            end

            rays_pos = rays(:,1:2); % Take only x and y position components
            rays_dir = rays(:,4:6);

            %% Rotate the generated points and vectors around the center point of the domain so that they are in the appropriate plane
            % https://math.stackexchange.com/questions/2093314/rotation-matrix-of-rotation-around-a-point-other-than-the-origin
            center_point = obj.size_VS/2; % Minimum of domain is always 0 in all 3 axes
            switch obj.boundary % Define rotation matrix based on boundary that flux is entering through
                case 'left' % Transform so that points lie in y-z plane, positive z-component of direction becomes positive x-component
                    new_order = [3,1,2];
                    flip_vec = [1,1,1];
                    shift_vec = [0,0,0];
                    
                case 'right' % Transform so that points lie in y-z plane, positive z-component of direction becomes negative x-component
                             % Shift to +x boundary
                    new_order = [3,1,2];
                    flip_vec = [-1,1,1];
                    shift_vec = [obj.size_VS(1),0,0];

                case 'bottom' % Transform so that points lie in x-z plane, positive z-component of direction becomes positive y-component 
                    new_order = [2,3,1];
                    flip_vec = [1,1,1];
                    shift_vec = [0,0,0];

                case 'top' % Transform so that points lie in x-z plane, positive z-component of direction becomes negative y-component 
                           % Shift to +y boundary
                    new_order = [2,3,1];
                    flip_vec = [1,-1,1];
                    shift_vec = [0,obj.size_VS(2),0];

                case 'back' % No transform required, this is the x-y plane with positive z-component pointing into the domain
                   new_order = [1,2,3];
                   flip_vec = [1,1,1];
                   shift_vec = [0,0,0];

                case 'front' % Flip z coordinate and shift to +z boundary
                   new_order = [1,2,3];
                   flip_vec = [1,1,-1];
                   shift_vec = [0,0,0];
            end
            rays_pos = [rays_pos,zeros(N_rays,1)];
            rays_pos = rays_pos(:,new_order);
            rays_pos = rays_pos+shift_vec;
            rays_dir = rays_dir(:,new_order).*flip_vec;
            
        end
    end
end

