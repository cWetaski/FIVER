function [VS_T_new, VS_dQ, VS_dT, count_itr] = equilibriumCondRad(N_rays,VS_T,VS_T_fixed,voxel_spaces,varargin)
% CONDRADEUILIBRIUM Solves the transient coupled heat transfer problem with conduction and radiation until equilibrium
% Inputs:
%   N_rays (1D double) [-]:                                 Number of rays per iteration, or an increasing vector of
%                                                           ray numbers. The solver will automatically go to the
%                                                           next ray level when it detects that it has reached a
%                                                           stationary solution. Unless N_rays is a scalar, the last
%                                                           level will only be iterated on 3 times. 
%   VS_T (3D double (T >= 0)) [K]:                          Initial temperature field of the voxel space
%   VS_T_fixed (3D logical):                                Logical array defining which voxels (if any) have fixed temperatures
%   voxel_spaces (1D cell of VoxelSpaces                    Can also be a singular VoxelSpace object.
%                                                           Each Voxel space object with the following:
%       opaque_voxels (3D logical):                             Stores which voxels are opaque (i.e., not PM, and not empty) voxels.
%       opaque_emissivities (3D double (0 <= eps <= 1)) [-]:    Stores the emissivity of opaque voxels
%       PM_absorption_coeffs (3D double (k >= 0)) [1/vx]:       Stores the linear absorption coefficient of PM voxels
%       surface_normals (3D cell):                              Stores 1x3 (normalized) surface normals for opaque surfaces
%                                                               Computed using GetNormalsAndSurfaceAreas.m
%       surface_areas (3D double (area >= 1) [vx^2]:            Stores area estimates for each opaque surface voxel
%                                                               Computed using GetNormalsAndSurfaceAreas.m
%       refractive_indexes (3D double (nn >= 1)) [-]:           Stores refractive indexes of medium
%       size (1x3 double (int) (sz >= 1)):                      Size of voxel space
%       voxel_scale (scalar double) [m/vx]:                     Scale of voxels
%       reflective_BCs (2x3 logical):                           Boundary conds: Rows are lower/upper bound, cols are XYZ. A
%                                                               reflective boundary reflects the ray specularly
%   k (scalar double) [W/(m-K)]:                            Thermal conductivity. Technically this value is the thermal
%                                                           thermal diffusivity (k/(rho*cp)) where rho = cp = 1, but
%                                                           it's an equilibrium solver so rho and cp don't affect the
%                                                           solution.
%   max_itr (scalar double (int)):                          Max number of iterations before returning.
%   spectral_bands (1D double (optional)) [um]:             Optional input for spectral bands. If spectral bands are
%                                                           included, then voxel_space must be a cell array of VoxelSpace 
%                                                           objects. Length of voxel_space array should be 1 longer than
%                                                           the length of spectral_bands variable (0 and inf are
%                                                           implied).
% Outputs:
%   VS_T_new (3D double):               Final temperature field
%   residuals_steady (3D double):       Residual of final temperature field from steady state heat transfer problem    

    %% Default Optional Parameters
    default_max_itrs = 100; % Maximum number of iterations per ray level
    default_final_level_itrs = 3; % Number of iterations at the final level (only used if N_rays has more than 1 value)
    default_N_prev = 5; % Number of previous iterations back to compare against for automatic stopping criteria
    default_C_converge = 1.1; % Scaling constant for convergence critera
    default_relax = 1; % Under/overrelaxation constant (i.e. <1 is underrelax and >1 is overrelax)
    default_spectral_bands = []; % Array of spectral band boundaries. 
    default_output_mode = 'concise';
    valid_output_modes = {'quiet','concise','verbose'};

    %% Parsing Varargin
    parser = inputParser;
 
    addParameter(parser,'MaxIterations',default_max_itrs);
    addParameter(parser,'NFinalLevelIterations',default_final_level_itrs);
    addParameter(parser,'NPreviousIterations',default_N_prev);
    addParameter(parser,'ConvergenceConstant',default_C_converge);
    addParameter(parser,'RelaxationConstant',default_relax);
    addParameter(parser,'SpectralBands',default_spectral_bands);
    addParameter(parser,'OutputMode',default_output_mode, ...
        @(x) any(validatestring(x,valid_output_modes)))

    parse(parser,varargin{:});
    
    max_itr = parser.Results.MaxIterations;
    N_prev = parser.Results.NPreviousIterations;
    final_level_itrs = parser.Results.NFinalLevelIterations;
    C_converge = parser.Results.ConvergenceConstant;
    C_relax = parser.Results.RelaxationConstant;
    spectral_bands = parser.Results.SpectralBands;
    output_mode = parser.Results.OutputMode;

    if strcmp(class(voxel_spaces),"VoxelSpace") 
        voxel_spaces = {voxel_spaces}; % allows this to work with a single voxel space as an input for gray radiation
    end

    %% Check for pure radiation
    if voxel_spaces{1}.thermal_conductivity == 0
        [VS_T_new,VS_dQ,VS_dT,count_itr] = equilibriumRad(N_rays,VS_T,VS_T_fixed,voxel_spaces,varargin{:}); % Go to pure radiation solver
        return
    end

    %% Constants
    %sigma = 5.670374419*10^(-8); % [W/(m^2-K^4)]: Stefan-Boltzmann constant
    
    %% Unpack relevant params from voxel space
    N_bands = length(spectral_bands)+1;
    size_VS = voxel_spaces{1}.size; 
    vx_scale = voxel_spaces{1}.voxel_scale;
    VS_alpha = voxel_spaces{1}.thermal_conductivity;
    max_PM_kappa = 0;
    external_flux_tot = 0;
    for i = 1:N_bands
        VS_PM_kappa = voxel_spaces{i}.PM_absorption_coeffs;
        max_PM_kappa = max([max_PM_kappa;VS_PM_kappa(:)]);
        for f = voxel_spaces{i}.external_fluxes
            external_flux_tot = external_flux_tot + f.power/(max(size_VS)*vx_scale).^2; % W/m2
        end
    end
    % if max_PM_kappa > 0
    %     N_param = max_PM_kappa/vx_scale*max(VS_alpha(:))/(4*sigma*max(VS_T(:))^3); % Estimate conduction-to-radiation factor (see Modest pp 725)
    %     dt = t_step_scaling*1./max_PM_kappa^2*N_param*vx_scale^2/(2*max(VS_alpha(:))); % [s]: Time step as N_param gets smaller, time step also gets smaller since more risk of diverging when radiation dominates, time step decreases with increasing PM_kappa since solution is exploding otherwise
    % else
    %     dt = t_step_scaling*(2*max(VS_alpha(:))/vx_scale)/(4*sigma*max(VS_T(:)).^3)*vx_scale^2; % No participating media
    % end
    %% Define conduction problem;
    thermal_mesh = createMesh3D(size_VS(1),size_VS(2),size_VS(3),size_VS(1)*vx_scale,size_VS(2)*vx_scale,size_VS(3)*vx_scale);
    BC = createBC(thermal_mesh); % all Neumann boundary condition structure is default
    
    alpha_cell = createCellVariable(thermal_mesh, VS_alpha); % assign the thermal diffusivity to the cells
    alpha_face = harmonicMean(alpha_cell); % calculate harmonic average of the diffusion coef on the cell faces
    alpha_face.xvalue(isnan(alpha_face.xvalue))=0;% Get rid of nans that arise from adjacent cells having 0 diffusivity!
    alpha_face.yvalue(isnan(alpha_face.yvalue))=0; 
    alpha_face.zvalue(isnan(alpha_face.zvalue))=0; 

    %% Temperature arrays
    VS_T_old = cell(N_prev,1); % Will keep track of N_prev previous iterations
    VS_T_prev = VS_T;
    VS_T_new = VS_T;
    VS_T_old{1} = VS_T;
    
  
    %% Define cell variables for FVTools solver
    T_prev = createCellVariable(thermal_mesh, VS_T_prev,BC); % initial values
    source_const_cell = createCellVariable(thermal_mesh,0);
    source_lin_cell = createCellVariable(thermal_mesh,0); % Initialize source term cell variable
    

    %% Set up fixed temperatures
    % Pad VS_T_spec so it is the same size as the cell variables (which have ghost cells);
    VS_T_spec_padded = padarray(VS_T_fixed,[1 1 1],0,'both');
    
    %% Get components of RHS and matrix of coefficients for diffusion and BCs (do not change on iterations)
    M_diff = diffusionTerm(alpha_face); % matrix of coefficients for the diffusion term
    [M_bc, RHS_bc] = boundaryCondition(BC); % matrix of coefficients and RHS vector for the BC
    M_size = size(M_diff);
    RHS_steady = RHS_bc;
    

    %% Steady state problem
    M_steady = M_bc - M_diff;
    
    %% Apply fixed coefficients
    fixed_inds = find(VS_T_spec_padded); % Get linear indices of fixed values
    M_fixed = sparse(fixed_inds,fixed_inds,1, M_size(1), M_size(2)); % Create a sparse matrix of 0s except for diagonals where temperature is fixed
    RHS_fixed = zeros(M_size(1),1);
    RHS_fixed(fixed_inds) = T_prev.value(fixed_inds); % Set fixed values on RHS
    
    % Replace entries with fixed values
    RHS_steady(fixed_inds) = 0;
    RHS_steady = RHS_steady + RHS_fixed;

    M_steady(fixed_inds,:) = 0;
    M_steady = M_steady + M_fixed;


    %% Initialize counters
    count_level = 1; % Counter which resets each time N_rays level is incremented
    count_itr = 0; % Total iteration counter
    count_final_level = 0; % Stop counter for when max level reached
    N_levels = length(N_rays); % Number of ray levels (assumed to be increasing)
    cur_level = 1; % Current ray level
    
    %% loop
    while count_final_level < final_level_itrs
        count_level = count_level + 1;
        if count_final_level > 0
            count_final_level = count_final_level + 1;
        end
        count_mod = mod(count_level-1,N_prev)+1; % Get iteration mod (1,2,...,N_prev,1,2,...,N_prev,1)
        [VS_dQ, VS_Q_emit_prev] = radiativeHeatFlowsMC(N_rays(cur_level),VS_T_prev,voxel_spaces,"SpectralBands",spectral_bands,"OutputMode",output_mode); % (3D Double) [W/vx^3]: Radiative flux divergence
        
        VS_dQ(VS_T_fixed) = 0;
        
        % Split dQ into a source term of form Y = aT + b; 
        source_const = (VS_dQ + 4*VS_Q_emit_prev)/vx_scale^3; % Divide by voxel element volume to get W/m3 i.e. source term comes from divergence of radiative heat flux
        source_const = padarray(source_const,[1,1,1],0,'both');
        source_const_cell.value = source_const;
        source_lin = 4*VS_Q_emit_prev./VS_T_prev/vx_scale^3;
        source_lin(isnan(source_lin)) = 0;
        source_lin = padarray(source_lin,[1,1,1],0,'both');
        source_lin_cell.value = source_lin; % Assign to source cell variable

        RHS_source_term = constantSourceTerm3D(source_const_cell); % Convert to RHS source term
        M_source_term = linearSourceTerm3D(source_lin_cell);

        M = M_steady+M_source_term;
        RHS = RHS_steady+RHS_source_term;
        
        % Set fixed temperature values;
        M(fixed_inds,:) = 0; % This should be redundant, haven't checked though if I can remove it safely yet
        M = M + M_fixed;

        RHS(fixed_inds) = 0;
        RHS = RHS + RHS_fixed;
        
        T_new = solvePDE(thermal_mesh,M, RHS); % Solving linear equation M*T=RHS

        VS_T_new = T_new.value(2:(end-1),2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)        
        VS_dT = VS_T_new - VS_T_prev; % Difference from previous temperature field
        
        if C_relax ~= 1 % Apply under/overrelaxation
            VS_T_new = VS_T_prev + C_relax*(VS_T_new-VS_T_prev);
        end
        VS_T_new(VS_T_new<0) = 0; % Do not allow negative temperatures

        % Adaptive stopping criteria
        % Premise:  Temperature changes at equilibrium will be only due to stochastic process, therefore residuals from
        %           the previous iteration will be of same order as residuals from many iterations ago

        if count_level > N_prev % The first N_prev iterations at a level, we just store the temperature fields
            
            R_sq_dT = sum(VS_dT.^2,'all'); % Sum of square residual from 0

            %VS_2_dT = VS_T_new - VS_T_old{mod(count_mod-3,N_prev)+1}; % Difference from 2 iterations ago
            VS_N_dT = VS_T_new - VS_T_old{count_mod}; % Temperature field difference from N_prev iterations ago
            R_sq_N_dT = sum(VS_N_dT.^2,'all'); % Sum of square residual from 0
            %R_sq_2_dT = sum(VS_2_dT.^2,'all'); 
            R_sq_ratio = R_sq_N_dT/R_sq_dT; % Ratio of square residuals.
            
            if strcmp(output_mode,'verbose')
                fprintf("Ratio of square residuals: %0.3f \n",R_sq_ratio)
            end

            if R_sq_ratio < C_converge || count_level > max_itr % stationary condition (slightly more permissive than pure radiation due to possible periodic behaviour
                cur_level = cur_level + 1;
                count_level = 0; % reset count (want at least N_prev iterations at each level).
                if cur_level == N_levels
                    count_final_level = 1; % Begin counting in last level until stopping criteria
                elseif cur_level > N_levels % Edge case for when there is only a single level
                    count_final_level = final_level_itrs; % Stop immediately (there is no higher level to proceed to)
                end
            end
        end
        % Update temperature fields
        VS_T_prev = VS_T_new;
        VS_T_old{count_mod} = VS_T_new;
    end
end

