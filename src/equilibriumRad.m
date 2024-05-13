function [VS_T_new, VS_dQ, VS_dT, count_itr, total_rays] = equilibriumRad(N_rays,VS_T,VS_T_fixed,voxel_space,varargin)
% MC_EQUILIBRIUM Iteratively finds the equilibrium temperature field for a voxel space
%   Conditions: 
%       Gray surfaces
%       Gray nonscattering medium with constant refractive index
%   Terms:
%       PM = Participating Media
%       VS = Voxel Space
%
% INPUTS:
%   N_rays (1D double) [-]:                                 Number of rays per iteration, or an increasing vector of
%                                                           ray numbers. The solver will automatically go to the
%                                                           next ray level when it detects that it has reached a
%                                                           stationary solution. Unless N_rays is a scalar, the last
%                                                           level will only be iterated on 3 times. 
%
%   VS_T (3D double (T >= 0)) [K]:                          Initial temperature field of voxel space
%
%   VS_T_fixed (3D logical):                                Logical array of voxels with fixed temperature (if any)
%                                                           The corresponding value is assumed to be the one given
%                                                           in the temperature_field of voxel_space.   
%
%   voxel_space (VoxelSpace object):                        Voxel space object (see src/VoxelSpace.m for fields)
%
%   varargin (optional inputs) enter as 'name', value pairs:
%       'MaxIterations' (scalar double (int)) (default = 200):  Maximum number of iterations per level.
%
%       'NFinalLevelIterations' (scalar double) (default = 3):  Number of iterations in the final level
%
%       'NPreviousIterations' (scalar double (int)) (default = 1):  Number of previous iterations stored for automatic stopping criteria 
%
%       'ConvergenceConstant' (scalar double) (default = 1):    Enter a value less than 1 for stricter convergence (but may
%                                                               never converge). Enter a value greater than 1 for easier
%                                                               convergence (but may stop iterating prematurely).
%
%       'RelaxationConstant' (scalar double) (default = 1):     Enter a value less than 1 for under relaxation, enter a
%                                                               a value greater than 1 for over relaxation
% 
%
%       'SpectralBandEdges' (1D double) [um] (default = []):    If spectral bands are included, then voxel_space must be a cell 
%                                                           array of VoxelSpace objects (each object provides the properties
%                                                           of the voxel space in its respective band). The length of 
%                                                           voxel_space array should be 1 longer than the length of 
%                                                           spectral_bands variable (0 and inf are implied).
% OUTPUTS:
%   VS_T (3D double) [K]:                                   Temperature field at equilibrium
%   VS_dQ (3D double) [W/m^3]:                              Radiative flux divergence field at equilibrium
%   VS_dT (3D double) [K]:                                  Change in temperature field at last iteration
%   count_itr (scalar):                                     Total number of iterations across all levels
%   total_rays (scalar):                                    Total number of rays traced

    %% Default Optional Parameters
    default_max_itrs = 100; % Maximum number of iterations per ray level
    default_final_level_itrs = 3; % Number of iterations at the final level (only used if N_rays has more than 1 value)
    default_N_prev = 5; % Number of previous iterations back to compare against for automatic stopping criteria
    default_C_converge = 1.05; % Scaling constant for convergence critera
    default_relax = 1; % Under/overrelaxation constant (i.e. <1 is underrelax and >1 is overrelax)
    default_spectral_band_edges = []; % Array of spectral band boundaries.
    default_output_mode = 'concise';
    valid_output_modes = {'quiet','concise','verbose'};
    
    %% Parsing Varargin
    parser = inputParser;
 
    addParameter(parser,'MaxIterations',default_max_itrs);
    addParameter(parser,'NFinalLevelIterations',default_final_level_itrs);
    addParameter(parser,'NPreviousIterations',default_N_prev);
    addParameter(parser,'ConvergenceConstant',default_C_converge);
    addParameter(parser,'RelaxationConstant',default_relax);
    addParameter(parser,'SpectralBandEdges',default_spectral_band_edges);
    addParameter(parser,'OutputMode',default_output_mode, ...
        @(x) any(validatestring(x,valid_output_modes)))

    parse(parser,varargin{:});
    
    max_itr = parser.Results.MaxIterations;
    N_prev = parser.Results.NPreviousIterations;
    final_level_itrs = parser.Results.NFinalLevelIterations;
    C_converge = parser.Results.ConvergenceConstant;
    C_relax = parser.Results.RelaxationConstant;
    spectral_bands = parser.Results.SpectralBandEdges;
    output_mode = parser.Results.OutputMode;
    
    %% Remaining Preamble

    if strcmp(class(voxel_space),"VoxelSpace") 
        voxel_space = {voxel_space}; % allows this to work with a single voxel space as an input for gray radiation
    end

    N_bands = length(spectral_bands)+1; 
    sigma = 5.670374419*10^(-8); % [W/(m^2-K^4)]: Stefan-Boltzmann constant
    vx_scale = voxel_space{1}.voxel_scale;

    if N_bands == 1
        % Used for time stepping in single band case
        VS_opaq_eps = voxel_space{1}.opaque_emissivities;
        VS_surf_areas = voxel_space{1}.surface_areas;
        VS_PM_kappa = voxel_space{1}.PM_absorption_coeffs;
        VS_nn = voxel_space{1}.refractive_indexes;
    end
    
    %% Initialize temperature arrays

    VS_T_old = cell(N_prev,1); % Will keep track of N_prev previous iterations of temperature field
    VS_T_old{1} = VS_T;
    VS_T_prev = VS_T;
    VS_T_new = VS_T;

    %% Initialize counters
    count_level = 1; % Counter which resets each time N_rays level is incremented (initially 1, since we have the initial temperature field)
    count_itr = 0; % Total iteration counter
    count_final_level = 0; % Stop counter for when final level reached
    N_levels = length(N_rays); % Number of ray levels (assumed to be increasing)
    total_rays = 0;
    cur_level = 1; % Current ray level

    while count_final_level <= final_level_itrs
        total_rays = total_rays + N_rays(cur_level); % total ray counter
        count_itr = count_itr + 1;
        count_level = count_level + 1;
        if count_final_level > 0
            count_final_level = count_final_level + 1;
        end
        count_mod = mod(count_level-1,N_prev)+1; % Get iteration multiple (1,2,...,N_prev,1,2,...,N_prev,1)
        [VS_dQ, VS_Q_emit_no_self,VS_Q_self_absorb] = radiativeHeatFlowsMC(N_rays(cur_level),VS_T_prev,voxel_space, ...
            "SpectralBandEdges",spectral_bands,"OutputMode",output_mode); % Ray trace
        
        VS_dQ(VS_T_fixed) = 0; % Set it so fixed temperatures don't change
        VS_Q_emit_prev = VS_Q_emit_no_self + VS_Q_self_absorb; % Add self-absorptions back to emissive power
        VS_Q_emit_new = VS_Q_emit_prev+C_relax*VS_dQ; % Get new emissions
        VS_Q_emit_negative = VS_Q_emit_new<0;
        Q_emit_negative = sum(VS_Q_emit_new(VS_Q_emit_negative)); % Check for negative values
        while Q_emit_negative<0
            VS_Q_emit_positive = (VS_Q_emit_new>0) & ~VS_T_fixed;
            N_positive = sum(VS_Q_emit_positive(:));
            dQ_negative = Q_emit_negative/N_positive;
            VS_Q_emit_new(VS_Q_emit_negative) = 0; % Remove negative values -> this adds energy to the system
            VS_Q_emit_new(VS_Q_emit_positive) = VS_Q_emit_new(VS_Q_emit_positive) + dQ_negative; % Offset added energy by removing energy from emitting voxels (distributed equally)
            VS_Q_emit_negative = VS_Q_emit_new<0;
            Q_emit_negative = sum(VS_Q_emit_new(VS_Q_emit_negative));
        end

        
        if N_bands == 1
            % Calculate new temperature field based on new emissions (different formula for surfaces vs PM voxels)
            VS_PM = logical(VS_PM_kappa); % Should make this part of VoxelSpace object, to be honest.
            VS_T_new(VS_PM) = (VS_Q_emit_new(VS_PM)./(vx_scale.^2*4.*VS_nn(VS_PM).^2*sigma.*VS_PM_kappa(VS_PM))).^(1/4);
            VS_T_new(VS_surf_areas>0) = (VS_Q_emit_new(VS_surf_areas>0)./(vx_scale.^2.*VS_nn(VS_surf_areas>0).^2.*sigma.*VS_surf_areas(VS_surf_areas>0).*VS_opaq_eps(VS_surf_areas>0))).^(1/4);        
        else % Can't easily invert Q_emit = f(T^4) since f is nonlinear function depending on wavelength bands
            VS_T_new =  VS_T_prev.*(VS_Q_emit_new./VS_Q_emit_prev).^(1/4); % estimate a linear proportionality to T^4 (i.e., AT^4 = Q -> T2/T1 = (Q2/Q1)^(1/4))
            nan_vals = find(isnan(VS_T_new)); % Occurs when VS_Q_emit_prev == 0 which can happen when number of rays is small
            for i = 1:length(nan_vals) % Properly invert temperature at the specific voxel (this is expensive, so we don't want to do it a lot)
                VS_T_new(nan_vals(i)) = fzero(@(T) TotalEmissivePowerSingle(spectral_bands,voxel_space,T,nan_vals(i))-VS_Q_emit_new(nan_vals(i)),[0,max(VS_T_new(:))]);
            end
        end

        % under and over relaxation
        %if C_relax ~= 1
        %    VS_T_new = VS_T_prev + C_relax*(VS_T_new-VS_T_prev);
        %end
        VS_dT = VS_T_new - VS_T_prev; % Difference from previous temperature field

        % Adaptive stopping criteria
        % Premise:  Temperature changes at equilibrium will be only due to stochastic process, therefore residuals from
        %           the previous iteration will be of same order as residuals from many iterations ago

        if count_level > N_prev % The first N_prev iterations at a level, we just store the temperature fields
            R_sq_dT = sum(VS_dT.^2,'all'); % Sum of square residual from 0

            VS_N_dT = VS_T_new - VS_T_old{count_mod}; % Temperature field difference from N_prev iterations ago
            R_sq_N_dT = sum(VS_N_dT.^2,'all'); % Sum of square residual from 0
            R_sq_ratio = R_sq_N_dT/R_sq_dT; % Ratio of square residuals.
            
            if strcmp(output_mode,'verbose')
                fprintf("Ratio of square residuals: %0.3f \n",R_sq_ratio)
            end
            
            if R_sq_ratio < C_converge || count_level > max_itr % stationary condition
                cur_level = cur_level + 1;
                count_level = 0; % reset count (want at least N_prev iterations at each level).
                if cur_level == N_levels
                    count_final_level = 1;
                elseif cur_level > N_levels % When there is only a single level -> end immediately
                    count_final_level = final_level_itrs+1; % End iteration immediately (there is no higher level to proceed to)
                end
            end
        end
        VS_T_old{count_mod} = VS_T_new; % Store temperature field for stopping check
        VS_T_prev = VS_T_new; % update temperature field for next monte carlo sim
    end
end