function VS_T_new = transientCondRad(N_rays,VS_T,VS_T_fixed,voxel_spaces,N_time_steps,time_step,use_internal_itr,spectral_bands)
% CONDRADTRANSIENT Solves the transient coupled heat transfer problem with conduction and radiation
% Inputs:
%   N_rays (scalar double (int)):       Number of rays to trace at each iteration
%   VS_T (3D double (T >= 0)) [K]:      Initial temperature field of the voxel space
%   VS_T_fixed (3D logical):            Logical array defining which voxels (if any) have fixed temperatures
%   voxel_space (VoxelSpace object):                            Voxel space object with the following:
%       opaque_voxels (3D logical):                             Stores which voxels are opaque (i.e., not PM, and not empty) voxels.
%       opaque_emissivities (3D double (0 <= eps <= 1)) [-]:    Stores the emissivity of opaque voxels
%       PM_absorption_coeffs (3D double (k >= 0)) [1/vx]:       Stores the linear absorption coefficient of PM voxels
%       surface_normals (3D cell):                              Stores 1x3 (normalized) surface normals for opaque surfaces
%                                                               Computed using GetNormalsAndSurfaceAreas.m
%       surface_areas (3D double (area >= 1) [vx^2]:            Stores area estimates for each opaque surface voxel
%                                                               Computed using GetNormalsAndSurfaceAreas.m
%       thermal_conds (3D double) [W/m-K]:                      Stores thermal conductivity of each voxel (only used in conduction problem)
%       densities (3D double) [kg/m^3]:                         Stores density of each voxel (only used in transient problem)
%       specific_heats (3D double) [J/kg-K]:                    Stores specific heat of each voxel (only used in transient problem)
%       refractive_index (scalar double (nn >= 1)) [-]:         Refractive index of medium (only homogenous mediums allowed for 
%                                                               now, and Snell's law not considered)
%       size (1x3 double (int) (sz >= 1)):                      Size of voxel space
%       voxel_scale (scalar double) [m/vx]:                     Scale of voxels
%       reflective_BCs (2x3 logical):                           Boundary conds: Rows are lower/upper bound, cols are XYZ. A
%                                                               reflective boundary reflects the ray specularly
%   N_t_steps (scalar double (int)):        Number of time steps to compute
%   t_step (scalar double) [s]:             Time step size
%   use_internal_itr (bool):                Whether to internally iterate on each step
%   spectral_bands (1D double (opt)) [um]:  Optional input for spectral bands. If spectral bands are
%                                           included, then voxel_space must be a cell array of VoxelSpace 
%                                           objects. Length of voxel_space array should be 1 longer than
%                                           the length of spectral_bands variable (0 and inf are implied).
%
% Outputs:
%   VS_T_new (3D double):               Final temperature field
%   residuals_steady (3D double):       Residual of final temperature field from steady state heat transfer problem    
    %% Check for spectral bands
    if nargin == 7
        spectral_bands = [];
    end
    if strcmp(class(voxel_spaces),"VoxelSpace") 
        voxel_spaces = {voxel_spaces}; % allows this to work with a single voxel space as an input for gray radiation
    end    

    %% Param
    internal_itr_num = 3; % Number of internal iterations if use_internal_itr is true
    % Can change this setting (eventually make it an argument with varargin along spectral bands!)


    % Unpack parameters for transient conduction problem take from the first voxel space since these properties can't
    % vary spectrally, therefore these fields only need to be defined in the first voxel space of the array for spectral
    % problem. 
    VS_kk = voxel_spaces{1}.thermal_conductivity; % [W/(m-K)]: Thermal conductivity
    VS_rho = voxel_spaces{1}.density; % [kg/m^3]: Density
    VS_cc = voxel_spaces{1}.specific_heat; % [J/(kg-K)]: Specific heat capacity
    
    VS_alpha = VS_kk./(VS_rho.*VS_cc); % [m^2/s]: Thermal diffusivity
    VS_alpha(isnan(VS_alpha)) = 0;
    size_VS = voxel_spaces{1}.size; 
    vx_scale = voxel_spaces{1}.voxel_scale;
    
    %% Define conduction problem;
    thermal_mesh = createMesh3D(size_VS(1),size_VS(2),size_VS(3),size_VS(1)*vx_scale,size_VS(2)*vx_scale,size_VS(3)*vx_scale);
    BC = createBC(thermal_mesh); % all Neumann boundary condition structure is default
    
    alpha_cell = createCellVariable(thermal_mesh, VS_alpha); % assign the thermal diffusivity to the cells
    alpha_face = harmonicMean(alpha_cell); % calculate harmonic average of the diffusion coef on the cell faces
    
    % Get components of RHS and matrix of coefficients for diffusion and BCs (do not change on iterations)
    M_diff = diffusionTerm(alpha_face); % matrix of coefficients for the diffusion term
    M_diff(isnan(M_diff)) = 0; % Nan Values to 0
    [M_bc, RHS_bc] = boundaryCondition(BC); % matrix of coefficients and RHS vector for the BC

    M_size = size(M_diff);
    
    %% define initial values
    VS_T_new = VS_T;
    VS_T_old = VS_T;
    T_old = createCellVariable(thermal_mesh, VS_T_new,BC); % initial values
    source_const_cell = createCellVariable(thermal_mesh,0); % Initialize source term cell variable
    source_lin_cell = createCellVariable(thermal_mesh,0);
    %% Set up fixed temperatures
    % Pad VS_T_spec so it is the same size as the cell variables (which have ghost cells);
    VS_T_spec_padded = padarray(VS_T_fixed,[1 1 1],0,'both');
    fixed_inds = find(VS_T_spec_padded); % Get linear indices of fixed values
    M_fixed = sparse(fixed_inds,fixed_inds,1, M_size(1), M_size(2)); % Create a sparse matrix of 0s except for diagonals where temperature is fixed
    RHS_fixed = zeros(M_size(1),1);
    RHS_fixed(fixed_inds) = T_old.value(fixed_inds); % Set fixed values on RHS
    
    M_steady = -M_diff+M_bc;
    M_steady(fixed_inds,:) = 0;
    M_steady = M_steady+M_fixed;
    RHS_steady = RHS_bc;
    RHS_steady(fixed_inds) = 0;
    RHS_steady = RHS_steady+RHS_fixed;


    %% loop
    for i = 1:N_time_steps
        internal_itr_count = 0;
        while internal_itr_count < internal_itr_num
            internal_itr_count = internal_itr_count + 1;
            heat_flows_tic = tic;
            [VS_dQ, VS_Q_emit_prev] = radiativeHeatFlowsMC(N_rays,VS_T_old,voxel_spaces,"SpectralBands",spectral_bands); % (3D Double) [W]: Power in/out of each voxel
            heat_flows_time = toc(heat_flows_tic);
            % Split dQ into a source term of form Y = aT + b; 
            source_const = (VS_dQ + 4*VS_Q_emit_prev)/(VS_cc.*VS_rho.*vx_scale^3); % K/s Divide by voxel element volume to get W/m3 i.e. source term comes from divergence of radiative heat flux
            source_const = padarray(source_const,[1,1,1],0,'both');
            source_const_cell.value = source_const;
            source_lin = 4*VS_Q_emit_prev./VS_T_prev/(VS_cc.*VS_rho.*vx_scale^3); % K/s
            source_lin(isnan(source_lin)) = 0;
            source_lin = padarray(source_lin,[1,1,1],0,'both');
            source_lin_cell.value = source_lin; % Assign to source cell variable
    
            RHS_source_term = constantSourceTerm3D(source_const_cell); % Convert to RHS source term
            M_source_term = linearSourceTerm3D(source_lin_cell);
    
            tic
            
            [M_trans, RHS_trans] = transientTerm(T_old, time_step); % Get transient term based on old temperature field. Don't quite understand how this works but it works
    
            M = M_steady + M_trans + M_source_term;
            RHS = RHS_steady + RHS_trans + RHS_source_term;
            
            % Set fixed temperature values
            M(fixed_inds,:) = 0;
            M = M + M_fixed;
    
            RHS(fixed_inds) = 0;
            RHS = RHS + RHS_fixed;

            T_new = solvePDE(thermal_mesh,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS
            VS_T_new = T_new.value(2:(end-1),2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)

            VS_T_old = VS_T_new; % Update temperature field for internal iteration (or time step, if internal iteration is finished)
            FVTool_time = toc;
            fprintf("Transient Internal Itr: heat flows time = %0.1f s    FVTool time = %0.1f s\n",heat_flows_time,FVTool_time)
            if ~use_internal_itr
                break
            end
        end
        % Update temperature field for FVTool transient term
        T_old = T_new;
    end
end

