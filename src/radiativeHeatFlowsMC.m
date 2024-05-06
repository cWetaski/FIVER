function [VS_Delta_Q, VS_Q_emit,VS_Q_absorb] = radiativeHeatFlowsMC(N_rays,VS_T,voxel_spaces,varargin)
% RADIATIVEHEATFLOWS Calculates the radiative heat flows of each voxel using monte carlo ray tracing.
%   N_rays are generated based on the temperature field and surface/PM properties. The rays are traced
%   until absorption or exiting the voxel space. All voxel spaces are assumed to have the same size. So many
%   arguments are used because MATLAB is slower at indexing large arrays when they are stored in structs/objects so
%   we're stuck with slightly overwhelming function definitions.
%   Conditions: 
%       Gray surfaces
%       Gray nonscattering medium with constant refractive index
%   Terms:
%       PM = Participating Media
%       VS = Voxel Space
%
% INPUTS:
%   N_rays (scalar double) [-]:                             Number of rays to generate
%   VS_T (3D double (T >= 0)) [K]:                          Temperature field of the voxel space
%   voxel_spaces (1D cell of VoxelSpace objects):           Voxel space object with the following:
%       opaque_voxels (3D logical):                             Stores which voxels are opaque (i.e., not PM, and not empty) voxels.
%       opaque_emissivities (3D double (0 <= eps <= 1)) [-]:    Stores the emissivity of opaque voxels
%       PM_absorption_coeffs (3D double (k >= 0)) [1/vx]:       Stores the linear absorption coefficient of PM voxels
%       surface_normals (3D cell):                              Stores 1x3 (normalized) surface normals for opaque surfaces
%                                                               Computed using GetNormalsAndSurfaceAreas.m
%       surface_areas (3D double (area >= 1) [vx^2]:            Stores area estimates for each opaque surface voxel
%                                                               Computed using GetNormalsAndSurfaceAreas.m
%       refractive_index (scalar double (nn >= 1)) [-]:         Refractive index of medium (only homogenous mediums allowed for 
%                                                               now, and Snell's law not considered)
%       size (1x3 double (int) (sz >= 1)):                      Size of voxel space
%       voxel_scale (scalar double) [m/vx]:                     (same for every voxel space) Scale of voxels 
%       reflective_BCs (2x3 logical):                           (same for every voxel space) Boundary conds: Rows 
%                                                               are lower/upper bound, cols are XYZ. Specular
%                                                               reflections
%  spectral_bands (1D double (optional)) [um]:              Optional input for spectral bands. If spectral bands are
%                                                           included, then voxel_space must be a cell array of VoxelSpace 
%                                                           objects. Length of voxel_space array should be 1 longer than
%                                                           the length of spectral_bands variable (0 and inf are
%                                                           implied).
% OUTPUTS:
%   Delta_Q (3D double) [W]:                                Net power in (positive) or out (negative) of each voxel
%   VS_Q_emit (3D double) [W]:                              Mathematical emissive power out of each voxel (used in
%                                                           equilibrium solver). Note this is not actually equal
%                                                           to the emissive power since actual ray emission locations 
%                                                           are generated probablistically.
%   VS_Q_absorb (3D double) [W]:                            Power absorbed by each voxel  
outer_timer = tic;
    %% Preamble
    
    default_spectral_bands = []; % Array of spectral band boundaries.
    default_output_mode = 'concise';
    valid_output_modes = {'quiet','concise','verbose'};

        %% Parsing Varargin
    parser = inputParser;
 
    addParameter(parser,'SpectralBands',default_spectral_bands);
    addParameter(parser,'OutputMode',default_output_mode, ...
        @(x) any(validatestring(x,valid_output_modes)))
    parse(parser,varargin{:});
    
    spectral_bands = parser.Results.SpectralBands;
    output_mode = parser.Results.OutputMode;

    if strcmp(class(voxel_spaces),"VoxelSpace") 
        voxel_spaces = {voxel_spaces}; % allows this to work with a single voxel space as an input for gray radiation
    end

    size_VS = voxel_spaces{1}.size;
    vx_scale = voxel_spaces{1}.voxel_scale; 
    reflective_BCs = voxel_spaces{1}.reflective_BCs; 
    
    N_vx_tot = prod(size_VS); % total number of voxels in domain
    sigma = 5.670374419*10^(-8); % [W/(m^2-K^4)]: Stefan-Boltzmann constant


    N_bands = length(spectral_bands)+1; % Get number of bands

    %% Initialize band-arrays
    VS_Q_emit = zeros(size_VS);
    VS_Q_emit_bands = cell(N_bands,1); 
    inds_vx_bands = cell(N_bands,1);
    Q_emit_vx_tot_bands = zeros(N_bands,1);
    N_vx_surf_bands = zeros(N_bands,1);
    N_vx_PM_bands = zeros(N_bands,1);
    Q_external_flux_band = zeros(N_bands,1);


    %% Get energy in each wavelength band
    for n = 1:N_bands
        % Unpack voxel space for band. Note that this does not pollute memory since nothing is modified. 

        VS_opaq_eps = voxel_spaces{n}.opaque_emissivities;
        VS_surf_areas = voxel_spaces{n}.surface_areas;
        VS_PM_kappa = voxel_spaces{n}.PM_absorption_coeffs;
        VS_nn = voxel_spaces{n}.refractive_indexes;
        external_fluxes = voxel_spaces{n}.external_fluxes;
 
        if N_bands == 1 % No need to compute spectral bands, just use n^2*sigma*T^4
            % Do not use logical indexing since in my tests it was slower to logically index 3-4 large matrices than
            % just multiplying them all as is.
            VS_Q_emit_surf_band = sigma*vx_scale^2*VS_surf_areas.*VS_opaq_eps.*VS_nn.*VS_T.^4; % (3D double) [W]: emissive power from each surface voxel within band
            VS_Q_emit_PM_band = 4*sigma*vx_scale^2*VS_PM_kappa.*VS_nn.*VS_T.^4; % (3D double) [W]: emissive power from each PM voxel within band (Modest eq 10.54)
        else    
            % Get lower and upper wavelength of each band
            if n == 1
                lb = 0;
                ub = spectral_bands(n);
            elseif n == N_bands
                lb = spectral_bands(n-1);
                ub = 10^8; 
            else
                lb = spectral_bands(n-1);
                ub = spectral_bands(n);
            end
            % Compute Powers (use logical indexing to reduce amount of computations for SpectralBandPower!)
            VS_Q_emit_surf_band = zeros(size_VS);
            VS_surfaces = logical(VS_surf_areas);
            VS_Q_emit_surf_band(VS_surfaces) = vx_scale.^2*VS_surf_areas(VS_surfaces).*VS_opaq_eps(VS_surfaces).*spectralBandPower(lb,ub,VS_T(VS_surfaces),VS_nn(VS_surfaces)); % (3D double) [W]: emissive power from each surface voxel within band according to Planck's law
            VS_Q_emit_PM_band = zeros(size_VS);
            VS_PM = logical(VS_PM_kappa);
            VS_Q_emit_PM_band(VS_PM) = 4*vx_scale.^2*VS_PM_kappa(VS_PM).*spectralBandPower(lb,ub,VS_T(VS_PM),VS_nn(VS_PM)); % (3D double) [W]: emissive power from each PM voxel within band (Modest eq 10.54)
        end

        Q_emit_surf_tot_band = sum(VS_Q_emit_surf_band,'all'); % (1D double) [W]: % Total emissive power from surface voxels
        Q_emit_PM_tot_band = sum(VS_Q_emit_PM_band,'all'); % (1D double) [W]: % Total emissive power from PM voxels
        Q_emit_vx_tot_bands(n) = Q_emit_surf_tot_band+Q_emit_PM_tot_band; % (1D double) [W]: Total emissive power in the band from voxels
        
        inds_vx_surf_band = find(VS_Q_emit_surf_band); % linear indices of emitting surface voxels in the band
        N_vx_surf_bands(n) = length(inds_vx_surf_band); % Number of emitting surface voxels in the band

        inds_vx_PM_band = find(VS_Q_emit_PM_band); % linear indices of emitting PM voxels in each band
        N_vx_PM_bands(n) = length(inds_vx_PM_band); % Number of emitting PM voxels in each band
        
        inds_vx_bands{n} = [inds_vx_surf_band;inds_vx_PM_band];
        VS_Q_emit_bands{n} = VS_Q_emit_surf_band+VS_Q_emit_PM_band;

        VS_Q_emit = VS_Q_emit + VS_Q_emit_bands{n};% Get total emissive power from each voxel
        
        for flux = external_fluxes
            Q_external_flux_band(n) = Q_external_flux_band(n)+flux.power; % Emissive power in band from external flux(es)
        end
    end



    %% Probablistically determine number of rays in each band
    Q_emit_tot_band = Q_emit_vx_tot_bands+Q_external_flux_band;
    Q_emit_tot = sum(Q_emit_tot_band); % (scalar double) [W]: Total emissive power across all bands
    % Get power power ray
    power_per_ray = Q_emit_tot/N_rays; % [W/ray]

    % Determine number of emissions from each wavelength band
    PDF_emit_band = (Q_emit_tot_band)/Q_emit_tot; % PDF of emissions by wavelength band
    CDF_emit_band = [0;cumsum(PDF_emit_band)]; % CDF of emissions by wavelength band
    
    N_rays_band = histcounts(rand(N_rays,1),CDF_emit_band)'; % Determine how many rays emitted from each band probablistically
    N_rays_vx_band = round(N_rays_band.*Q_emit_vx_tot_bands./(Q_emit_tot_band)); % Not probabilistic, number of rays from voxels in each band
    N_rays_external_flux_band = N_rays_band-N_rays_vx_band; % Number of rays from external flux in each band
    
    
    %% Initialize results arrays
    emission_counts = zeros(N_vx_tot,1); % Combined emission counts vector
    absorption_pos = cell(N_bands,1); % Absorptions positions

    %% Monitoring
    ray_tracing_time = zeros(N_bands,1);
    ray_generation_time = zeros(N_bands,1);

    for n = 1:N_bands % Don't run this in parallel, since some bands might have very little radiation thus those workers just wait around (or may have fewer bands than workers)
        if N_rays_band(n) == 0
            continue;
        end
        % Unpack parts of voxel space which are used for ray tracing for current band 
        VS_opaq = voxel_spaces{n}.opaque_voxels;
        VS_opaq_eps = voxel_spaces{n}.opaque_emissivities;
        VS_surf_norms = voxel_spaces{n}.surface_normals;
        %surf_norm_inds = voxel_spaces{n}.surface_normal_inds;
        %surf_norms_lin = voxel_spaces{n}.surface_normals_lin;
        VS_PM_kappa = voxel_spaces{n}.PM_absorption_coeffs;
        VS_nn = voxel_spaces{n}.refractive_indexes;
        external_fluxes = voxel_spaces{n}.external_fluxes;

        %% Handle all Emissions
        if N_rays_vx_band(n) > 0
            inds_vx_band = inds_vx_bands{n};  
            Q_emit_filter = VS_Q_emit_bands{n}(inds_vx_band); % Just take emitting voxels (becomes a 1D array)
            emission_counts_band = floor(N_rays_vx_band(n)*Q_emit_filter/Q_emit_vx_tot_bands(n)); % Assign rays deterministically
            N_rays_remaining = N_rays_vx_band(n)-sum(emission_counts_band);
            if N_rays_remaining > 0
                remaining_Q_emit = Q_emit_filter-emission_counts_band*power_per_ray; % Get remaining emissive power
                if any(remaining_Q_emit(:)<0)
                    remaining_Q_emit(remaining_Q_emit<0) = 0;
                end
                CDF_emit_vx = [0;cumsum(remaining_Q_emit/(sum(remaining_Q_emit)))]; % CDF of remaining emissions from voxels in band
                
                emission_counts_band = emission_counts_band + histcounts(rand(N_rays_remaining,1),CDF_emit_vx)'; % Assign remaining emissions probablistically
            end
            % Filter 0 emissions voxels
            N_vx_surf_band = nnz(emission_counts_band(1:N_vx_surf_bands(n))); % nnz -> count nonzero elements
            N_vx_PM_band = nnz(emission_counts_band((N_vx_surf_bands(n)+1):end));
            N_vx_band = N_vx_surf_band + N_vx_PM_band;
            inds_vx_band(emission_counts_band == 0) = [];
            emission_counts_band(emission_counts_band==0) = [];
           
            emission_counts(inds_vx_band) = emission_counts(inds_vx_band) + emission_counts_band; % Increment total emissions 
            
            surf_norms_band = zeros(N_vx_band,3);% Even though only the first N_vx_surf_band elements are used, this array must be the size of N_vx_band or else parfor throws an error
            surf_norms_band(1:N_vx_surf_band,:) = cell2mat(VS_surf_norms(inds_vx_band(1:N_vx_surf_band))); % Get surface normals for emitting surface voxels
             [X, Y, Z] = ind2sub(size_VS,inds_vx_band); % Get subscript positions of emitters
            
            tic
            rays_pos = cell(N_vx_band,1);
            rays_dir = cell(N_vx_band,1);        
            for i = 1:N_vx_band % I used to parfor this but for some reason it randomly (like 1 in 3000 times) throws an feval error. It also doesn't save much time anyway and is in fact slower in some cases.
                cur_N_rays = emission_counts_band(i);
                cur_vx = [X(i),Y(i),Z(i)];
                if i <= N_vx_surf_band % Surface voxel
                    % Generate rays
                    [rays_pos{i}, rays_dir{i}] = generateSurfaceRays(cur_vx,surf_norms_band(i,:),cur_N_rays);
                else % PM voxel
                    % generate N_rays random positions within voxel
                    rays_pos{i} = [cur_vx(1),cur_vx(2),cur_vx(3)] - rand(cur_N_rays,3);
                    
                    % Generate N_emit random directions
                    phis = 2*pi*rand(cur_N_rays,1); % Modest eq 21.13b
                    thetas = acos(1-2*rand(cur_N_rays,1)); % Modest eq 21.13c
                    rays_dir{i} = [sin(thetas).*cos(phis),sin(thetas).*sin(phis),cos(thetas)]; % Modest eq 21.14
                end
            end
            rays_pos = cell2mat(rays_pos);
            rays_dir = cell2mat(rays_dir);
        else
            rays_pos = [];
            rays_dir = [];
        end
        if N_rays_external_flux_band > 0
            N_external_flux = length(external_fluxes);
            N_rays_flux = zeros(N_external_flux,1);
            for ii = 1:N_external_flux
                if ii < N_external_flux
                    N_rays_flux(ii) = round(N_rays_external_flux_band(n)*external_fluxes(ii)/Q_external_flux_band(n));
                else
                    N_rays_flux(ii) = N_rays_external_flux_band(n)-sum(N_rays_flux);
                end
                [rays_pos_external_flux,rays_dir_external_flux] = external_fluxes(ii).GenerateRays(N_rays_flux(ii));
            end
            rays_pos = [rays_pos;rays_pos_external_flux]; %#ok<AGROW>
            rays_dir = [rays_dir;rays_dir_external_flux]; %#ok<AGROW>
        end
        %if any(size_VS==1) % This makes traversal faster should we have a 2D or 1D domain since we don't have to constantly bounce off the boundary
        %    rays_dir(:,size_VS==1) = 0;
        %    rays_dir = rays_dir./vecnorm(rays_dir,2,2);
        %end

        ray_generation_time(n) = toc;
        tic
        [absorption_pos_band,events] = traverseRays(rays_pos,rays_dir,VS_opaq,VS_opaq_eps,VS_surf_norms,VS_PM_kappa,VS_nn,reflective_BCs,size_VS);
        ray_tracing_time(n) = toc;

        absorption_pos_band(events == 1,:) = []; % Remove exit events
        absorption_pos{n} = absorption_pos_band;
        % Update absorption counts


    end % Iterating over spectral bands
    absorption_pos = cell2mat(absorption_pos);
    absorption_counts = histcounts(absorption_pos(:,1) + ...
                                   (absorption_pos(:,2) - 1)*size_VS(1) + ...
                                   (absorption_pos(:,3) - 1)*size_VS(1)*size_VS(2), ... # This is faster version of sub2ind
                                    1:(N_vx_tot+1))';
    %% Calculate Delta_Q
    % Net number of absorbed rays for each voxel * power per ray = Power out of each voxel
    VS_Q_absorb = absorption_counts*power_per_ray;
    VS_Q_absorb = reshape(VS_Q_absorb,size_VS);

    VS_Delta_Q = (absorption_counts - emission_counts)*power_per_ray; % (1D double) [W]
    VS_Delta_Q = reshape(VS_Delta_Q,size_VS); % (3D double) [W]: reshape to 3D voxel space
    
    total_time = toc(outer_timer);
    if any(contains({'concise','verbose'},output_mode))
        fprintf("Heat Flows: total time = %0.1f s    ray generation = %0.1f s   ray tracing = %0.1f s\n", ... 
                total_time,sum(ray_generation_time),sum(ray_tracing_time));  
    if true %% DEBUG
        VS_dQ_nan = VS_Delta_Q;
        VS_dQ_nan(VS_surf_areas==0)=nan;
        mean_dQ_z = squeeze(sum(VS_dQ_nan,[1,2],'omitnan'));

    end
end