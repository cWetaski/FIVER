classdef TES_System < matlab.mixin.Copyable
    %TES_SYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        voxel_spaces; % 1D cell array of Voxel Spaces
        spectral_band_edges; % 1D vector (um) defining band edges
        Rosseland_kappas; % 1D cell array of 3D arrays of linear absorption coefficients which are being Rosseland'd NOT IMPLEMENTED YET (POSSIBLY EVER)
        pipe_systems; % 1D cell array of Pipes
        heaters; % 1D cell of Heater object(s)
        mesh; % FVTool mesh
        VS_T; % Temperature field -> modified by charging/discharging
        matrix_ind_map; % 3D array mapping i,j,k coordinates (including ghost cells) to matrix rows/cols
        
        %%% TES GEOMETRY
        VS_TES; % 3D Logical array defining geometry of the TES
        VS_pipes; % 3D Logical array defining which voxels/grid cells are occupied by pipes through the TES. VS_TES is not set to false where VS_pipes is set to true
        VS_heaters; % 3D Logical array defining which voxels/grid cells are occupied by heaters in the TES. VS_TES is not set to false where VS_heaters is set to true
        VS_TES_material; % 3D Logical array defining the space actually occupied by TES material (i.e., excluding pipes and heaters)

        %%% PARTS OF SYSTEM OF EQUATIONS WHICH DON'T CHANGE BASED ON TEMPERATURE
        M_diff; % Diffusion matrix
        M_BC; % Boundary condition matrix
        RHS_BC; % Boundary condition RHS
        
        %%% PROPERTIES FOR FIXED TEMPERATURE
        VS_T_fixed; % 3D Logical array defining which temperature values should be held fixed
        fixed_indices; % (1D vector) corresponding to rows in the matrices which have fixed temperatures
        M_fixed; % Fixed temperatures matrix
        RHS_fixed; % Fixed temperatures RHS
        
        %%% PHYSICAL PROPERTIES OF THE MEDIUM
        VS_alpha; % 3D array of thermal diffusivity
        TES_material TES_Material; % TES_Material object
        
        %%% METADATA
        size_VS; % (1x3 double (int)) Size of the voxel space (copied from voxel_spaces)
        size_M; % (scalar) Size of the matrices (=== prod(size_VS+2))
        vx_scale; % (1x3 double) m/voxel length of each voxel
        energy_capacity;
        minimum_temperature = 30+273.15; % 30°C -> can set to a different value, but then you should call updateEnergyCapacity()

        %%% OTHER
        figure_handle; % Handle to figure
        last_discharge_power; %
    end
    
    methods
        function obj = TES_System(VS_TES,material_name,vx_scale,T)
            %TES_SYSTEM Construct an instance of this class
            %   Detailed explanation goes here  
            obj.VS_TES = VS_TES;
            obj.TES_material = TES_Material(material_name);

            [Nx,Ny,Nz] = size(VS_TES);
            obj.size_VS = [Nx,Ny,Nz];            
            obj.size_M = (Nx+2)*(Ny+2)*(Nz+2);

            obj.mesh = createMesh3D(Nx,Ny,Nz,Nx*vx_scale(1),Ny*vx_scale(2),Nz*vx_scale(3));
            obj.matrix_ind_map = reshape(1:obj.size_M, Nx+2, Ny+2,Nz+2);

            alpha = obj.TES_material.thermal_conductivity/(obj.TES_material.specific_heat*obj.TES_material.density);
            obj.VS_alpha = zeros(Nx,Ny,Nz);
            obj.VS_alpha(VS_TES) = alpha;
            obj.addDiffusionMatrix();
            
            obj.vx_scale = vx_scale;
            obj.VS_T_fixed = false(Nx,Ny,Nz);
            obj.VS_pipes = false(Nx,Ny,Nz);
            obj.VS_heaters = false(Nx,Ny,Nz);
            obj.VS_TES_material = VS_TES;

            obj.updateEnergyCapacity();

            if nargin == 2 % T is optional argument
                obj.VS_T = zeros(Nx,Ny,Nz);
            else
                obj.VS_T = ones(Nx,Ny,Nz).*T; % T can be scalar or 3D array
            end
        end

        function addVoxelSpaces(obj,voxel_spaces,spectral_band_edges)
            if ~isa(voxel_spaces,'cell')
                error("Voxel Spaces must be cell array of voxel spaces (even if non-spectral, in which case there should be a singule cell element)")
            end
            if nargin == 3
                obj.spectral_band_edges = spectral_band_edges;
            end
            obj.voxel_spaces = voxel_spaces;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% addFixedTemperature

        function addFixedTemperature(obj,VS_T_fixed,VS_T) % VS
            obj.VS_T_fixed = VS_T_fixed;
            fixed_inds = find(padarray(VS_T_fixed,[1 1 1],0,'both')); % Get linear indices of fixed values, padding for ghost cells
            obj.fixed_indices = fixed_inds;
            obj.M_fixed = sparse(fixed_inds,fixed_inds,1, obj.size_M, obj.size_M); % Create a sparse matrix of 0s except for diagonals where temperature is fixed
            
            obj.RHS_fixed = zeros(obj.size_M,1);
            if nargin == 2 % VS_T is optional, if not included, uses the existing temperature field
                VS_T_padded = padarray(obj.VS_T,[1 1 1],0,'both');
            else
                VS_T_padded = padarray(VS_T,[1 1 1],0,'both');
            end
            obj.RHS_fixed(fixed_inds) = VS_T_padded(fixed_inds); % Set fixed values on RHS
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% addDiffusionMatrix

        function addDiffusionMatrix(obj,VS_alpha)
            if nargin == 2                
                if size(VS_alpha) ~= obj.size_VS
                    error('VS_alpha does not match size of domain')
                end
                obj.VS_alpha = VS_alpha;
            end

            alpha_cell = createCellVariable(obj.mesh, obj.VS_alpha); % assign the thermal diffusivity to the cells
            alpha_face = harmonicMean(alpha_cell); % calculate harmonic average of the diffusion coef on the cell faces
            alpha_face.xvalue(isnan(alpha_face.xvalue))=0;% Get rid of nans that arise from adjacent cells having 0 diffusivity!
            alpha_face.yvalue(isnan(alpha_face.yvalue))=0; 
            alpha_face.zvalue(isnan(alpha_face.zvalue))=0;
            obj.M_diff = diffusionTerm(alpha_face); % matrix of coefficients for the diffusion term
            
            for pipe_system = obj.pipe_systems % If pipe systems already added to system, this deletes the appropriate rows from the diffusion matrix
                pipe_inds = vertcat(pipe_system.inds_pipe_boundary,pipe_system.inds_pipe_interior);
                for nn = 1:size(pipe_inds,1)
                    iii = squeeze(obj.matrix_ind_map(inds(nn,1),inds(nn,2),2:end-1));
                    obj.M_diff(iii,:) = 0; % Delete diffusion term from pipe locations
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% addDefaultBC

        function addDefaultBC(obj)
            BC = createBC(obj.mesh); % all Neumann (==0) boundary condition structure is default
            [obj.M_BC, obj.RHS_BC] = boundaryCondition(BC); % matrix of coefficients and RHS vector for the BC
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% addPipe

        function addPipe(obj,x_loc,y_loc,shape,dim_spec,roughness,fluid_property_table)
            switch shape
                case "circle"
                    radius_vx = dim_spec;
                    pipe_system = CirclePipe(x_loc,y_loc,radius_vx,roughness,obj.vx_scale,obj.size_VS,fluid_property_table);  
                case "square"
                    xy_width_vx = dim_spec;
                    pipe_system = SquarePipe(x_loc,y_loc,xy_width_vx,roughness,obj.vx_scale,obj.size_VS,fluid_property_table);
            end
            pipe_inds = vertcat(pipe_system.inds_pipe_boundary,pipe_system.inds_pipe_interior);
            for nn = 1:size(pipe_inds,1)
                iii = squeeze(obj.matrix_ind_map(pipe_inds(nn,1)+1,pipe_inds(nn,2)+1,2:end-1)); % add 1 due to ghost cells
                obj.M_diff(iii,:) = 0; % Delete diffusion term from pipe locations
            end
            obj.pipe_systems{end+1} = pipe_system;
            VS_pipe_3D = repmat(pipe_system.VS_pipe_2D,[1,1,obj.size_VS(3)]);
            obj.VS_pipes(VS_pipe_3D) = true;
            obj.VS_TES_material(VS_pipe_3D) = false;
            obj.removeFromVoxelSpaces(VS_pipe_3D);
            obj.updateEnergyCapacity();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% addHeater

        function addHeater(obj,x_loc,y_loc,z_bottom,z_top,radius,void_radius,emissivity,max_temperature,fixed_T_value)
            if nargin == 8
                fixed_T_value = 0; % fixed_T_value is an optional input for what temperature should be outputed in the heater's location
                % In the future, it's possible we solve directly for the heater's temperature.
            end
            obj.heaters{end+1} = Heater(x_loc,y_loc,z_bottom,z_top,radius,void_radius,emissivity,max_temperature,obj.vx_scale,obj.size_VS);
            if true % If heater size is significant (in voxel units)
                %% Add a void in TES where the heater is
                VS_heater = false(obj.size_VS);
                for x_dif = -ceil(void_radius):ceil(void_radius)
                    if x_dif+x_loc == obj.size_VS(1) % Edge case exactly on boundary
                        x_ind = obj.size_VS(1);
                    else
                        x_ind = floor(x_dif+x_loc)+1; %
                    end
                    if x_ind < 1 || x_ind > obj.size_VS(1)
                        continue
                    end
                    cur_x_center = x_ind-0.5; % Center x-location of current voxel
                    for y_dif = -ceil(void_radius):ceil(void_radius)
                        if y_dif+y_loc == obj.size_VS(2)
                            y_ind = obj.size_VS(2);
                        else
                            y_ind = floor(y_dif+y_loc)+1;
                        end
                        if y_ind < 1 || y_ind > obj.size_VS(2)
                            continue
                        end
                       cur_y_center = y_ind-0.5; % Center y-location of current voxel
                       cur_center_distance = ((cur_x_center-x_loc)^2+(cur_y_center-y_loc)^2)^(0.5); % Distance from center of heater to center of current voxel
                       cur_radial_vec = [cur_x_center,cur_y_center]-[x_loc,y_loc];
                       cur_radial_vec = cur_radial_vec/norm(cur_radial_vec);
                     
                       if cur_center_distance <= void_radius+0.5/max(abs(cur_radial_vec))
                           VS_heater(x_ind,y_ind,:) = true;
                       end
                    end
                end
                obj.VS_heaters(VS_heater) = true;
                obj.VS_alpha(VS_heater) = 0;
                obj.VS_T_fixed(VS_heater) = true;
                obj.VS_T(VS_heater) = fixed_T_value;
                obj.VS_TES_material(VS_heater) = false;
                obj.removeFromVoxelSpaces(VS_heater);
                obj.addDiffusionMatrix(); % Recompute the diffusion matrix
                obj.addFixedTemperature(obj.VS_T_fixed);
                obj.updateEnergyCapacity()
            end % Otherwise, we assume the heater is infinitesimal and has no effect on the overall TES properties
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getHeatCapacity
        
        function C = getHeatCapacity(obj) % big C heat capacity
            C = sum(obj.VS_TES_material*prod(obj.vx_scale)*obj.TES_material.density*obj.TES_material.specific_heat,'all'); % J/K
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% updateEnergyCapacity
        
        function updateEnergyCapacity(obj) % big C heat capacity
            obj.energy_capacity = obj.getHeatCapacity()*(obj.TES_material.annealing_temperature-obj.minimum_temperature);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% createVoxelSpaces

        function createVoxelSpaces(obj,spectral_band_edges,ns,specular_BCs)
            % Note, this function lacks customizability (e.g. with respect to emissivity, multiple materials, refractive interfaces etc.), can expand/modify in the future
            N_band = max(length(spectral_band_edges)-1,1);
            obj.spectral_band_edges = spectral_band_edges;

            VS_opaq_PM_band = ~obj.VS_TES;
            VS_opaq_opaq_band = ~obj.VS_pipes & ~obj.VS_heaters;
            all_zeros = zeros(obj.size_VS);
            all_ones = ones(obj.size_VS);
            obj.voxel_spaces = cell(N_band,1);            
            [VS_surf_norms_PM,VS_surf_areas_PM] = getNormalsAndSurfaceAreas(VS_opaq_PM_band,obj.vx_scale,ns);
            [VS_surf_norms_opaq,VS_surf_areas_opaq] = getNormalsAndSurfaceAreas(VS_opaq_opaq_band,obj.vx_scale,ns);
            for nn = 1:N_band
                voxel_space = VoxelSpace();
                if isempty(spectral_band_edges)
                    band_l = -1; 
                    band_u = -1;
                else
                    band_l = spectral_band_edges(nn);
                    band_u = spectral_band_edges(nn+1);
                end
                PM_kappa = obj.TES_material.getAbsorptivity(band_l,band_u); % 1/m
                if PM_kappa ~= -1
                    voxel_space.opaque_voxels = VS_opaq_PM_band;
                    voxel_space.opaque_emissivities = all_zeros;
                    voxel_space.surface_areas = VS_surf_areas_PM;
                    voxel_space.surface_normals = VS_surf_norms_PM;
                    VS_PM_kappa = all_zeros;
                    %PM_kappa = 1e9;
                    VS_PM_kappa(obj.VS_TES_material) = PM_kappa; % 1/vx absorption coefficients
                    voxel_space.PM_absorption_coeffs = VS_PM_kappa;
                    n_refrac = obj.TES_material.getRefractiveIndex(band_l,band_u);
                    voxel_space.refractive_indexes = all_ones*n_refrac;
                    if ~isempty(spectral_band_edges) %PM_kappa*vx_scale > 2 % (CURRENTLY NO ROSSELAND) Optical thickness of 1 voxel
                        T_range = linspace(300,obj.TES_material.annealing_temperature,100);
                        k_Rosseland_range = 4/3/PM_kappa*spectralBandDerivative(band_l,band_u,T_range,n_refrac); % Nominal emissive power in the band
                        fprintf("Max k_Rosseland [%0.2f,%0.2f] = %0.4f W/m-K, (kappa = %0.1f 1/m)\n",band_l,band_u,max(k_Rosseland_range(:)),PM_kappa)
                        
                        if max(k_Rosseland_range(:)) < 0.1*obj.TES_material.thermal_conductivity && PM_kappa*obj.vx_scale(1)>1e-1
                            PM_kappa = -1;
                            fprintf("Setting this band to opaque")
                        end
                        %k_Rosseland = k_Rosseland + k_Rosseland_cur; % Add to Rosseland conductivity
                        % PM_kappa == -1; % Currently no Rosseland
                        %fprintf("vx_kappa: %0.2f, nom_emissive_power: %0.2f, k_Rosseland_cur: %0.2f \n",PM_kappa*vx_scale,nom_emissive_power,k_Rosseland_cur)
                    end
                end
                if PM_kappa == -1
                    voxel_space.opaque_voxels = VS_opaq_opaq_band;
                    voxel_space.opaque_emissivities = double(VS_opaq_opaq_band);
                    voxel_space.opaque_emissivities(obj.VS_TES_material) = obj.TES_material.getEmissivity(obj.VS_T(obj.VS_TES_material));
                    voxel_space.PM_absorption_coeffs = all_zeros;
                    voxel_space.surface_areas = VS_surf_areas_opaq;
                    voxel_space.surface_normals = VS_surf_norms_opaq;
                    voxel_space.refractive_indexes = all_ones;
                end
                voxel_space.size = obj.size_VS;
                voxel_space.voxel_scale = obj.vx_scale;
                voxel_space.ns_normals = ns;
                voxel_space.thermal_conductivity = obj.TES_material.thermal_conductivity;
                voxel_space.density = obj.TES_material.density;
                voxel_space.specific_heat = obj.TES_material.specific_heat;
                voxel_space.reflective_BCs = specular_BCs;
                obj.voxel_spaces{nn} = voxel_space;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% chargeAndDischarge

        function output = chargeAndDischarge(obj,time_step_size,N_time_steps,N_rays,T_fluid_i,T_process,nominal_discharge_power,nominal_heater_power_fraction)
            T_old = createCellVariable(obj.mesh,obj.VS_T);
            VS_T_new = obj.VS_T;
            source_const_cell = createCellVariable(obj.mesh,0);
            source_lin_cell = createCellVariable(obj.mesh,0);
            N_charge = length(obj.heaters);
            
            Q_rad_t = zeros(N_time_steps,1);
            Q_flow_t = zeros(N_time_steps,1);
            p_loss_t = zeros(N_time_steps,1);
            m_dot_t = zeros(N_time_steps,1);
            T_fluid_o_t = zeros(N_time_steps,1);
           
            heat_capacity = obj.TES_material.density*obj.TES_material.specific_heat;   
            vol_heat_capacity = heat_capacity*prod(obj.vx_scale); % [J/K] % heat capacity of single volume element
                        
            dT_max = 0;
            N_pipe = length(obj.pipe_systems);


            for tt = 1:N_time_steps
                VS_T_old = VS_T_new;
                VS_T_cur = VS_T_new;
                T_max = max(VS_T_old(obj.VS_TES_material));
                if nominal_heater_power_fraction > 0
                    if T_max > obj.TES_material.annealing_temperature
                        disp("TES temperature exceeds material annealing temperature, heater power set to 0")
                        heater_power_fraction = 0;
                    else
                        heater_power_fraction = nominal_heater_power_fraction;
                    end
                else
                    heater_power_fraction = 0;
                end

                [M_trans, RHS_trans] = transientTerm(T_old, time_step_size); % Get transient term based on old temperature field

                q_rel_diff = 1e6;
                q_rel_diff2 = 1e6;
                q_rel_diff_change = -1e6;
                count = 0;

                while q_rel_diff > 0.005 && q_rel_diff2 > 1e-8 && count < 10 && abs(q_rel_diff_change) > 0.001*q_rel_diff
                    count = count+1;
                    q_rel_diff_prev = q_rel_diff;
                    %% Radiation
                    if N_rays > 0
                        for nn = 1:N_charge
                            if heater_power_fraction > 0
                                obj.heaters{nn}.turnOnHeater(obj,heater_power_fraction)
                            else
                                obj.heaters{nn}.last_power_fraction = 0;
                            end
                        end
                        N_rays_rem = N_rays;
                        VS_dQ = 0; VS_Q_emit_no_self_abs = 0;
                        while N_rays_rem > 0
                            if N_rays_rem >= 20e6
                                N_rays_cur = 20e6;
                            else
                                N_rays_cur = N_rays_rem;
                            end
                            [VS_dQ_cur,VS_Q_emit_no_self_abs_cur_rays] = radiativeHeatFlowsMC(N_rays_cur,VS_T_cur,obj.voxel_spaces,'SpectralBandEdges',obj.spectral_band_edges);
                            VS_dQ = VS_dQ + VS_dQ_cur*N_rays_cur/N_rays;
                            VS_Q_emit_no_self_abs = VS_Q_emit_no_self_abs + VS_Q_emit_no_self_abs_cur_rays*N_rays_cur/N_rays;
                            N_rays_rem = N_rays_rem-N_rays_cur;
                        end
                            
                        if heater_power_fraction>0
                            for nn = 1:N_charge
                                obj.heaters{nn}.turnOffHeater(obj)
                            end
                        end
                        VS_dQ(obj.VS_T_fixed) = 0;

                        Q_rad_t(tt) = sum(VS_dQ(:));
                        for nn = 1:N_pipe
                            obj.pipe_systems{nn}.S_wall_z = obj.pipe_systems{nn}.getWallVariable(VS_dQ/vol_heat_capacity);
                        end
                        
                        if true
                            % Split dQ into a source term of form Y = aT + b;
                            % We don't use the "real linear approximation" because that has very slow convergence, owing to the fact
                            % that it only accounts for emissions from a voxel,
                            lin_fac = min(1,50/max(obj.VS_T(:))); lin_fac = 0;
                            source_const = (VS_dQ + 4*VS_Q_emit_no_self_abs.*lin_fac)/vol_heat_capacity; % units of Kelvin/s -> constant heat generation term  
                            %source_const = VS_dQ/vol_heat_capacity;
                            source_const(isnan(source_const)) = 0;
                            source_const = padarray(source_const,[1,1,1],0,'both');
                            source_const_cell.value = source_const;
                            source_lin = 4*lin_fac.*VS_Q_emit_no_self_abs./VS_T_cur/vol_heat_capacity; % units of 1/s -> linear heat generation term
                            source_lin(isnan(source_lin)) = 0;
                            %source_lin(:) = 0;
                            source_lin = padarray(source_lin,[1,1,1],0,'both');
                            source_lin_cell.value = source_lin; % Assign to source cell variable
                        else
                            source_const = VS_dQ/vol_heat_capacity;
                            source_const = padarray(source_const,[1,1,1],0,'both');
                            source_const_cell.value = source_const;
                            source_lin_cell.value = zeros(obj.size_VS+2);
                        end                
                    else
                        for nn = 1:N_pipe
                            obj.pipe_systems{nn}.S_wall_z = zeros(obj.size_VS(3),1);
                        end
                    end

                    
                        
                    RHS_source_term = constantSourceTerm3D(source_const_cell); % Convert to RHS source term
                    M_source_term = linearSourceTerm3D(source_lin_cell);

                RHS_pipe = zeros(obj.size_M,1);
                M_pipe = sparse(obj.size_M,obj.size_M);
                discharge_power = nominal_discharge_power;
                if N_pipe > 0
                    if nominal_discharge_power > 0
                        % The discharge power can "wiggle" in transition regime, therefore, fsolve can get stuck, so we implement brent's method directly.
                        % It is easy to find a bounding interval since discharge power goes to 0 as m_dot goes to 0 and gets very large as m_dot gets very large
                        [~,cp_in] = obj.pipe_systems{1}.getFluidProperties(T_fluid_i);
                        [~,cp_out_T_process] = obj.pipe_systems{1}.getFluidProperties(T_process);
                        [~,cp_out_T_max] = obj.pipe_systems{1}.getFluidProperties(T_max);

                        if count > 1
                            p_adjust = p_adjust + (abs(Q_TES_t(tt))-abs(Q_net(tt)))*0.5;
                        else
                            p_adjust = 0;
                        end
                        m_dot_max = (nominal_discharge_power-p_adjust)/(cp_out_T_process*T_process-cp_in*T_fluid_i);
                        m_dot_min = (nominal_discharge_power-p_adjust)/(cp_out_T_max*T_max-cp_in*T_fluid_i);
                        
                        error_fun = @(m_dot) obj.getDischargePower(m_dot,T_fluid_i)-nominal_discharge_power+p_adjust;
                        error_max = error_fun(m_dot_max);
                        if error_max < 0
                            m_dot_cur = m_dot_max;
                            discharge_power = nominal_discharge_power+error_max;
                        else
                            tol = 1e-6*nominal_discharge_power; % Tolerance defined relative to desired discharge power
                            m_dot_cur = fsolveBrent(error_fun,m_dot_min,m_dot_max,tol);
                        end
                    else
                        m_dot_cur = 0;
                    end
                end
                Q_flow = zeros(N_pipe,1);
                T_fluid_o = zeros(N_pipe,1);
                p_loss = zeros(N_pipe,1);
                for nn = 1:N_pipe
                    pipe_system = obj.pipe_systems{nn};
                    inds = vertcat(pipe_system.inds_pipe_interior,pipe_system.inds_pipe_boundary)+1; % plus 1 to account for ghost cells
                    for mm = 1:size(inds,1)
                        iii = squeeze(obj.matrix_ind_map(inds(mm,1),inds(mm,2),2:end-1));
                        M_trans(sub2ind(size(M_trans),iii,iii)) = 0; % Delete transient term from pipe locations
                        RHS_trans(iii) = 0;
                    end
                    [h_mean_z,T_fluid_z,Q_flow(nn),T_fluid_o(nn),p_loss(nn)] = solveFlow(obj,pipe_system,m_dot_cur,T_fluid_i);

                    %fprintf("Number of iterations for heat transfer coefficient: %d\n",count);
                    T_wall_z = pipe_system.getWallTemperature(obj.VS_T,obj.TES_material.thermal_conductivity,h_mean_z,T_fluid_z);
                    fprintf("CHECK Q: %0.4f kW (should be %0.4f kW) \n", 4*pipe_system.A_cs/pipe_system.D_hydraulic*(sum(h_mean_z.*(T_wall_z-T_fluid_z))/1000)*obj.vx_scale(3),nominal_discharge_power/1000);
                    [M_pipe_cur, RHS_pipe_cur] = pipe_system.getPipeEquation(obj.mesh,obj.TES_material.thermal_conductivity,h_mean_z,T_fluid_z);    
                    M_pipe = M_pipe+M_pipe_cur;
                    RHS_pipe = RHS_pipe + RHS_pipe_cur;
                end
                air_table = readtable('air_properties.txt','Delimiter',' '); % Table of air properties vs T at p=1atm              
                p_loss_t(tt) = sum(p_loss);
                T_fluid_o_t(tt) = mean(T_fluid_o);
                Q_flow_t(tt) = -sum(Q_flow);

                    %% Define system of equations
                    M = M_trans - obj.M_diff + M_pipe + obj.M_BC + M_source_term;
                    RHS = RHS_trans + RHS_pipe + obj.RHS_BC + RHS_source_term;
                    
                    % Set fixed temperature values;
                    M(obj.fixed_indices,:) = 0;
                    M = M + obj.M_fixed;
                
                    RHS(obj.fixed_indices) = 0;
                    RHS = RHS + obj.RHS_fixed;

                    Q_net(tt) = Q_rad_t(tt)-discharge_power;

                
                    T_new = solvePDE(obj.mesh,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS but returns a cell variable
                    VS_T_new = T_new.value(2:(end-1),2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)
                    VS_T_cur = VS_T_cur + (VS_T_new-VS_T_cur)*1; % underrelax internal iteration!
                    %VS_T_plot = VS_T_new;
                    %VS_T_plot(~obj.VS_TES_material) = nan;
                    %T_slice = VS_T_plot(:,:,1);
                    %h = imagesc(T_slice-273.15);
                    %set(h, 'AlphaData', ~isnan(T_slice))
                    %colorbar();
                    %clim([1100,1300])
                    %obj.VS_T = VS_T_cur;
                    Delta_E_t(tt) = sum((VS_T_new-VS_T_old).*(obj.VS_TES_material)*vol_heat_capacity,'all'); % Joules
                    Q_TES_t(tt) = Delta_E_t(tt)/time_step_size;
                   
                    q_diff = abs(Q_TES_t(tt)-Q_net(tt));
                    q_rel_diff = q_diff/abs(Q_net(tt));
                    fprintf("q_rel_diff = %0.4f, q_diff = %0.4f kW\n",q_rel_diff,(abs(Q_TES_t(tt))-abs(Q_net(tt)))/1000)
                    q_rel_diff2 = (q_diff)/obj.energy_capacity;
                    q_rel_diff_change = (q_rel_diff-q_rel_diff_prev);
                    if any(VS_T_new(:)<0)
                       debug = 0
                    end
                end
                obj.VS_T = VS_T_new; % Update temperature
                    
                m_dot_t(tt) = m_dot_cur;

                if q_rel_diff > 0.001
                    %fprintf("q_rel_diff = %0.3f \n",q_rel_diff)
                end
                cur_T_max = max(VS_T_new(obj.VS_TES_material));
                cur_T_min = min(VS_T_new(obj.VS_TES_material));
                T_range = [cur_T_min,cur_T_max];
                cur_dT_max = cur_T_max-cur_T_min;
                if cur_dT_max > dT_max
                    dT_max = cur_dT_max;
                end
            end
            output = struct();
            output.TES_energy_change = Delta_E_t;
            output.TES_power_rate = Q_TES_t;
            output.T_fluid_out = T_fluid_o_t-273.15; % in Celsius
            output.mass_flow_rate_per_pipe = m_dot_t; % Mass flow rate IN EACH PIPE, not total!
            output.pressure_loss = p_loss_t;
            output.TES_dT_max = dT_max;
            output.TES_T_range = T_range-273.15; % in Celsius
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getDischargePower

        function Q_tot = getDischargePower(obj,m_dot,T_fluid_i)
            N_pipe = length(obj.pipe_systems);
            Q_tot = 0;
            for nn = 1:N_pipe
                [~,~,Q_cur,~] = solveFlow(obj,obj.pipe_systems{nn},m_dot,T_fluid_i);
                Q_tot = Q_tot + Q_cur;
            end
        end

        function [h_mean_z,T_fluid_z,Q_flow,T_fluid_o,p_loss] = solveFlow(obj,pipe_system,m_dot,T_fluid_i)
            T_wall_z = pipe_system.last_T_wall_z;
            if isempty(pipe_system.last_T_wall_z)
                T_wall_z = ones(obj.size_VS(3),1)*mean(obj.VS_T(obj.VS_TES & ~obj.VS_pipes & ~obj.VS_heaters));
            end
            max_change = 1e6;
            count = 0;
            while(max_change>0.01) && count < 10
                count = count + 1;
                T_wall_z_old = T_wall_z;
               
                [h_mean_z,T_fluid_z,Q_flow,T_fluid_o,p_loss] = pipe_system.getHeatTransferCoefficient(T_wall_z,T_fluid_i,m_dot);
                T_wall_z = pipe_system.getWallTemperature(obj.VS_T,obj.TES_material.thermal_conductivity,h_mean_z,T_fluid_z);
                max_change = max(abs(T_wall_z - T_wall_z_old));       
            end  
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% visualizeTES

        function fig_handle = visualizeTES(obj,z_coord)
            if ~isempty(obj.figure_handle)
                set(0, 'current', obj.figure_handle);
                fig_handle = obj.figure_handle;
                hold off;
            else
                fig_handle = figure;
                obj.figure_handle = fig_handle;
            end
            T_slice = obj.VS_T(:,:,z_coord);
            T_slice(~obj.VS_TES_material(:,:,z_coord)) = nan;
            [max_domain_length,ind_max] = max([obj.size_VS(1)*obj.vx_scale(1),obj.size_VS(2)*obj.vx_scale(2)]);
            
            if max_domain_length > 1
                units = 'm';
                multiplier = 1*obj.vx_scale(ind_max);
            elseif max_domain_length > 0.01
                units = 'cm';
                multiplier = 100*obj.vx_scale(ind_max);
            else
                units = 'mm';
                multiplier = 1000*obj.vx_scale(ind_max);
            end
            himage = imagesc(T_slice-273.15,'XData',[0.5,obj.size_VS(1)-0.5]*multiplier,'YData',[0.5,obj.size_VS(2)-0.5]*multiplier);
            set(himage, 'AlphaData', ~isnan(T_slice))
            xlabel(sprintf("x (%s)",units));
            ylabel(sprintf("y (%s)",units));
            set(gca, 'YDir','normal');
            cb = colorbar();
            ylabel(cb,'T (°C)','FontSize',16,'Rotation',0)
           
            T_min = floor((min(obj.VS_T(obj.VS_TES_material))-273.15)/100)*100;
            T_max = ceil((max(obj.VS_T(obj.VS_TES_material))-273.15)/100)*100;
            clim([T_min,T_max])
            
            hold on
            
            N_heaters = length(obj.heaters);
            loc_heaters = zeros(N_heaters,2);
            total_heater_power = 0;
            for nn = 1:N_heaters
                loc_heaters(nn,:) = [obj.heaters{nn}.x_loc*multiplier,obj.heaters{nn}.y_loc*multiplier];
                total_heater_power = total_heater_power + obj.heaters{nn}.max_power*obj.heaters{nn}.last_power_fraction;
                heater_temp = obj.heaters{nn}.max_temperature*(obj.heaters{nn}.last_power_fraction)^(1/4);
            end
            h_heaters = scatter(loc_heaters(:,1),loc_heaters(:,2),10,'red');
            
            legend(h_heaters,sprintf("Heaters: \n Q_{tot} = %0.0f W \n T = %0.0f °C)",total_heater_power,heater_temp-273.15));
            drawnow();
            hold off
        end

        function removeFromVoxelSpaces(obj,VS)
            N_band = length(obj.voxel_spaces);
            for nn = 1:N_band
                if any(obj.voxel_spaces{nn}.PM_absorption_coeffs(VS))
                    obj.voxel_spaces{nn}.PM_absorption_coeffs(VS) = 0; % potentially assign absortion based on fluid properties
                end
                if any(obj.voxel_spaces{nn}.opaque_voxels(VS))
                    obj.voxel_spaces{nn}.opaque_voxels(VS) = 0;
                    obj.voxel_spaces{nn}.opaque_emissivities(VS) = 0;
                    [obj.voxel_spaces{nn}.surface_normalsVS_surf_norms,obj.voxel_spaces{nn}.surface_areas] = getNormalsAndSurfaceAreas(obj.voxel_spaces{nn}.opaque_voxels,obj.voxel_spaces{nn}.ns_normals);
                end
            end
        end
        
        function applyRosseland(obj,optical_thickness_threshold) % DO NOT USE
            N_band = length(obj.voxel_spaces);
            if isempty(obj.Rosseland_kappas)
                obj.Rosseland_kappas = cell(N_band,1);
            end
            for nn = 1:N_band
                VS_Rosseland_kappa = obj.voxel_spaces{nn}.PM_absorption_coeffs;
                VS_threshold = VS_Rosseland_kappa*obj.vx_scale>optical_thickness_threshold;
                if ~any(VS_threshold)
                    continue
                end
                VS_Rosseland_kappa(~VS_threshold) = 0;
                obj.Rosseland_kappas{nn} = VS_Rosseland_kappa;
                obj.voxel_spaces{nn}.PM_absorption_coeffs(VS_threshold) = 0;
                obj.voxel_spaces{nn}.opaque_voxels(VS_threshold) = true;
                obj.voxel_spaces{nn}.opaque_emissivities(VS_threshold) = 1;
                [obj.voxel_spaces{nn}.VS_surf_norms,obj.voxel_spaces{nn}.VS_surf_areas] = getNormalsAndSurfaceAreas(obj.voxel_spaces{nn}.opaque_voxels,obj.voxel_spaces{nn}.ns_normals);
            end
        end
        function removeRosseland(obj) % DO NOT USE
            N_band = length(obj.voxel_spaces);
            if isempty(obj.Rosseland_kappas)
                return
            end
            for nn = 1:N_band
                obj.voxel_spaces{nn}.PM_absorption_coeffs(obj.Rosseland_kappas>0) = obj.Rosseland_kappas(obj.Rosseland_kappas>0);
                obj.voxel_spaces{nn}.opaque_voxels(obj.Rosseland_kappas>0) = false;
                obj.voxel_spaces{nn}.opaque_emissivities(obj.Rosseland_kappas>0) = 0;
                [obj.voxel_spaces{nn}.VS_surf_norms,obj.voxel_spaces{nn}.VS_surf_areas] = getNormalsAndSurfaceAreas(obj.voxel_spaces{nn}.opaque_voxels,obj.voxel_spaces{nn}.ns_normals);
            end  
            obj.Rosseland_kappas = []; % Clear VS_Rosselands
        end
        function VS_Rosseland_k = getRosselandConductivity(obj) % DO NOT USE
            if isempty(obj.Rosseland_kappas)
                VS_Rosseland_k = 0;
            else
                VS_Rosseland_k = zeros(obj.size_VS);
                N_band = length(obj.Rosseland_kappas);
                for nn = 1:N_band
                    band_l = obj.spectral_band_edges(nn);
                    band_u = obj.spectral_band_edges(nn+1);
                    VS_Rosseland = logical(obj.Rosseland_kappas{nn});
                    VS_black_body_cur = spectralBandPower(band_l,band_u,obj.VS_T(VS_Rosseland),obj.voxel_spaces{nn}.refractive_indexes(VS_Rosseland));
                    VS_Rosseland_k = VS_Rosseland_k + 16/3*VS_black_body_cur./obj.VS_T./obj.Rosseland_kappas{nn};
                end
            end
        end
    end
end

