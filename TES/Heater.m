classdef Heater < matlab.mixin.Copyable
    %HEATER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius; % Radius in voxel units
        void_radius; % Radius of the void in the TES (must be greater than the radius)
        x_loc; % x location in voxel coordinates
        y_loc; % y location in voxel coordinates
        z_bottom; % z location of bottom of heater
        z_top; % z location of top of heater
        surface_area; % Surface area and power together determine the temperature (assuming blackbody) and therefore spectral distribution
        emissivity; % Blackbody emissivity
        max_temperature; % Maximum surface temperature
        max_power; % Maximum output power (assuming no incoming radiation)
        ray_gen_function;     
        last_power_fraction; % For plotting, the power fraction of the heater's last use.
    end
    
    methods
        function obj = Heater(x_loc,y_loc,z_bottom,z_top,radius,void_radius,emissivity,max_temperature,vx_scale,size_VS)
            %HEATER Construct an instance of this class
            %   Detailed explanation goes here
            L = (z_top - z_bottom)*vx_scale(3); % [m]
            obj.x_loc = x_loc;
            obj.y_loc = y_loc;
            obj.z_bottom = z_bottom;
            obj.z_top = z_top;
            obj.radius = radius;
            obj.void_radius = void_radius;
            obj.surface_area = 2*radius*vx_scale(1)*pi*L; % [m2] Only consider curved surface of cylinder, assumes x and y voxel_scale are the same!!
            obj.emissivity = emissivity;
            obj.max_temperature = max_temperature;
            obj.max_power = obj.surface_area*emissivity*5.67e-8*max_temperature^4;

            %% Check boundary
            % If it's located on the boundary, we assume we only want half the radiation
            % Therefore, we bump the heater off the boundary a very small amount, and halve the power.
            % The distribution remains the same but we obtain the correct amount of radiation as long as the corresponding boundary has
            % a reflective boundary condition.
            if x_loc == 0
                x_loc = x_loc + 1e-10; % take off boundary
                obj.max_power = obj.max_power/2; % divide power by 2 since half the power should be leaving the edge
                obj.surface_area = obj.surface_area/2;
            elseif x_loc == size_VS(1)
                x_loc = x_loc - 1e-10;
                obj.max_power = obj.max_power/2; % divide po
                obj.surface_area = obj.surface_area/2;
            end
            if y_loc == 0
                y_loc = y_loc + 1e-10; % take off boundary
                obj.max_power = obj.max_power/2; %
                obj.surface_area = obj.surface_area/2;
            elseif y_loc == size_VS(2)
                y_loc = y_loc - 1e-10;
                obj.max_power = obj.max_power/2; %
                obj.surface_area = obj.surface_area/2;
            end

            %% Get ray function
            % The code below is written for a general line flux between two points. It could
            % be simplified quite a bit given that the current Heater definition only allows axis parallel to z-axis
            % Not sure if it has significant affect on performance.
            p1 = [x_loc,y_loc,z_bottom];
            p2 = [x_loc,y_loc,z_top];
            line = (p2-p1);
            line_norm = line/norm(line); % Normalized
            line_orthog = null(line_norm)';
            line_orthog_vec = line_orthog(1,:);
            ray_pos_generator = @(N_rays) p1 + line.*rand(N_rays,1);
            lambertian_ray_dir = @(phis,thetas) [sin(thetas).*cos(phis), sin(thetas).*sin(phis), cos(thetas)];
            lambertian_ray_generator = @(N_rays) lambertian_ray_dir(2*pi*rand(N_rays,1),asin(sqrt(rand(N_rays,1))));
            rotate_about_u = @(rays,thetas,u) [rays(:,1).*(cos(thetas)+u(1)^2*(1-cos(thetas)))+rays(:,2).*(u(1)*u(2)*(1-cos(thetas))-u(3)*sin(thetas))+rays(:,3).*(u(1)*u(2)*(1-cos(thetas))+u(2)*sin(thetas)), ...
                                               rays(:,1).*(u(1)*u(2)*(1-cos(thetas))+u(3)*sin(thetas))+rays(:,2).*(cos(thetas)+u(2)^2*(1-cos(thetas)))+rays(:,3).*(u(2)*u(3)*(1-cos(thetas))-u(1)*sin(thetas)), ...
                                               rays(:,1).*(u(3)*u(1)*(1-cos(thetas))-u(2)*sin(thetas))+rays(:,2).*(u(3)*u(2)*(1-cos(thetas))+u(1)*sin(thetas))+rays(:,3).*(cos(thetas)+u(3)^2*(1-cos(thetas)))];
            transform_to_line = @(rays) transformCoord(rays,line_orthog_vec);
            ray_dir_generator = @(N_rays) rotate_about_u(transform_to_line(lambertian_ray_generator(N_rays)),2*pi*rand(N_rays,1),line_norm);
            
            
            obj.ray_gen_function = @(N_rays) [ray_pos_generator(N_rays),ray_dir_generator(N_rays)];
        end
        
        function turnOnHeater(obj,TES_system,power_fraction)
            %METHOD1 Summary of this method goes here
            %   Adds charging fluxes to the voxel spaces of the TES_system
            
            if isempty(TES_system.spectral_band_edges) % Gray radiation system
                TES_system.voxel_spaces{1}.fluxes{end+1} = ChargingFlux(obj.max_power*power_fraction,obj.ray_gen_function,obj);
            else
                N_band = length(TES_system.spectral_band_edges)-1;
                T = obj.max_temperature*(power_fraction^(1/4));
                for n = 1:N_band
                    band_l = TES_system.spectral_band_edges(n);
                    band_u = TES_system.spectral_band_edges(n+1);
                    band_power = spectralBandPower(band_l,band_u,T,1)*obj.surface_area*obj.emissivity; % W: Need to figure out refractive index here...
                    TES_system.voxel_spaces{n}.addFlux(ChargingFlux(band_power,obj.ray_gen_function,obj));
                end
            end
            obj.last_power_fraction = power_fraction;
        end

        function turnOffHeater(obj,TES_system)
            % Removes all ChargingFluxes associated with the heater from the target cell array of voxel_spaces
            N_band = length(TES_system.voxel_spaces);
            for n = 1:N_band
                cur_voxel_space = TES_system.voxel_spaces{n}; % note that voxel_spaces are handle objects, so this doesn't create a copy
                N_fluxes = length(cur_voxel_space.fluxes);
                N_delete = 0;
                for n2 = 1:N_fluxes
                    m = n2 - N_delete;
                    cur_flux = cur_voxel_space.fluxes{m};
                    if isa(cur_flux,"ChargingFlux") && cur_flux.heater == obj
                        cur_voxel_space.fluxes(m) = [];
                        N_delete = N_delete + 1;
                    end
                end
            end
        end
    end
end

