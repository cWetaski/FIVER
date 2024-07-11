%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07 (Charles Wetaski)

function [VS_dQ, VS_Q_emit_no_self,VS_Q_self_absorb] = radiativeHeatFlowsMC(N_rays,VS_T,voxel_spaces,varargin)
    % RADIATIVEHEATFLOWS Calculates the radiative heat flows of each voxel using monte carlo ray tracing.
    %   N_rays are generated based on the temperature field and voxel_spaces properties. The rays are traced
    %   until absorption or exiting the voxel space.
    %
    % INPUTS:
    %   N_rays (scalar double) [-]:                 Number of rays to generate
    %   VS_T (3D double (T >= 0)) [K]:              Temperature field of the voxel space
    %   voxel_spaces (1D cell of VoxelSpaces):      Can also be a singular VoxelSpace object. See VoxelSpace.m for properties
    %   varargin:                                   Optional ("Name",value) pair arguments, default values defined below.
    %       "SpectralBandEdges" (1D double) [um]:   List of wavelengths corresponding to spectral bands. Should have length
    %                                                   which is 1 greater than the length of the vector of voxel_spaces,
    %                                                   as the properties in each voxel_space correspond to the 
    %                                                   properties in the respective wavelength band. 
    %       "OutputMode"                            Determines how much info is written to the command window.
    %                                                   Options are ["quiet","concise","verbose"]      
    % OUTPUTS:
    %   VS_dQ (3D double) [W]:                      Net change in power of each voxel
    %   VS_Q_emit_no_self (3D double) [W]:          Emissive power of each voxel excluding self absorptions.
    %                                                   Self absorptions refer to rays which are emitted by
    %                                                   the same ray they are absorbed by.
    %   VS_Q_self_absorb (3D double) [W]:           Self-absorption power of each voxel.
    %
    outer_timer = tic;
    %% Preamble
    
    default_spectral_band_edges = []; % Array of spectral band boundaries.
    default_output_mode = 'concise';
    valid_output_modes = {'quiet','concise','verbose'};

        %% Parsing Varargin
    parser = inputParser;
 
    addParameter(parser,'SpectralBandEdges',default_spectral_band_edges);
    addParameter(parser,'OutputMode',default_output_mode, ...
        @(x) any(validatestring(x,valid_output_modes)))
    parse(parser,varargin{:});
    
    spectral_band_edges = parser.Results.SpectralBandEdges;
    output_mode = parser.Results.OutputMode;

    if strcmp(class(voxel_spaces),"VoxelSpace") 
        voxel_spaces = {voxel_spaces}; % allows this to work with a single voxel space as an input for gray radiation
    end
    
    Vxyz = voxel_spaces{1}.Vxyz;
    size_VS = voxel_spaces{1}.size;
    vx_scale = voxel_spaces{1}.voxel_scale; 
    reflective_BCs = voxel_spaces{1}.reflective_BCs; 
    
    N_vx_tot = prod(size_VS); % total number of voxels in domain
    sigma = 5.670374419*10^(-8); % [W/(m^2-K^4)]: Stefan-Boltzmann constant


    N_bands = max(1,length(spectral_band_edges)-1); % Get number of bands

    %% Initialize band-arrays
    VS_Q_emit = zeros(size_VS);
    VS_Q_emit_bands = cell(N_bands,1); 
    inds_vx_bands = cell(N_bands,1);
    Q_emit_vx_tot_bands = zeros(N_bands,1);
    N_vx_surf_bands = zeros(N_bands,1);
    N_vx_PM_bands = zeros(N_bands,1);
    Q_flux_band = zeros(N_bands,1);

    %% Get energy in each wavelength band
    for n = 1:N_bands
        % Unpack voxel space for band. Note that this does not pollute memory since nothing is modified. 

        VS_opaq_eps = voxel_spaces{n}.opaque_emissivities;
        VS_surf_areas = voxel_spaces{n}.surface_areas;
        VS_PM_kappa = voxel_spaces{n}.PM_absorption_coeffs;
        VS_nn = voxel_spaces{n}.refractive_indexes;
        fluxes = voxel_spaces{n}.fluxes;

        if isempty(spectral_band_edges) % No need to compute spectral bands, just use n^2*sigma*T^4
            % Do not use logical indexing since in my tests it was slower to logically index 3-4 large matrices than
            % just multiplying them all as is.
            VS_Q_emit_surf_band = sigma*vx_scale^2*VS_surf_areas.*VS_opaq_eps.*VS_nn.*VS_T.^4; % (3D double) [W]: emissive power from each surface voxel within band
            VS_Q_emit_PM_band = 4*sigma*prod(Vxyz)*vx_scale^2*VS_PM_kappa.*VS_nn.*VS_T.^4; % (3D double) [W]: emissive power from each PM voxel within band (Modest eq 10.54)
        else    
            % Get lower and upper wavelength of each band
            lb = spectral_band_edges(n);
            ub = spectral_band_edges(n+1);
            % Compute Powers (use logical indexing to reduce amount of computations for SpectralBandPower!)
            VS_Q_emit_surf_band = zeros(size_VS);
            VS_surfaces = logical(VS_surf_areas);
            VS_Q_emit_surf_band(VS_surfaces) = vx_scale.^2*VS_surf_areas(VS_surfaces).*VS_opaq_eps(VS_surfaces).*spectralBandPower(lb,ub,VS_T(VS_surfaces),VS_nn(VS_surfaces)); % (3D double) [W]: emissive power from each surface voxel within band according to Planck's law
            VS_Q_emit_PM_band = zeros(size_VS);
            VS_PM = logical(VS_PM_kappa);
            VS_Q_emit_PM_band(VS_PM) = 4*prod(Vxyz)*vx_scale^2*VS_PM_kappa(VS_PM).*spectralBandPower(lb,ub,VS_T(VS_PM),VS_nn(VS_PM)); % (3D double) [W]: emissive power from each PM voxel within band (Modest eq 10.54)
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
        
        for m = 1:length(fluxes)
            Q_flux_band(n) = Q_flux_band(n)+fluxes{m}.power; % Emissive power in band from flux(es)
        end
    end



    %% Probablistically determine number of rays in each band
    Q_emit_tot_band = Q_emit_vx_tot_bands+Q_flux_band;
    Q_emit_tot = sum(Q_emit_tot_band); % (scalar double) [W]: Total emissive power across all bands
    % Get power power ray
    if Q_emit_tot == 0
        VS_dQ = zeros(size_VS);
        VS_Q_emit_no_self = zeros(size_VS);
        VS_Q_self_absorb = zeros(size_VS);
        return
    end
    power_per_ray = Q_emit_tot/N_rays; % [W/ray]

    % Determine number of emissions from each wavelength band
    PDF_emit_band = (Q_emit_tot_band)/Q_emit_tot; % PDF of emissions by wavelength band
    CDF_emit_band = [0;cumsum(PDF_emit_band)]; % CDF of emissions by wavelength band
    
    %N_rays_band = histcounts(rand(N_rays,1),CDF_emit_band)'; % Determine how many rays emitted from each band probablistically
    N_rays_band = floor(N_rays*(Q_emit_tot_band)/Q_emit_tot);
    N_rays_remaining = N_rays-sum(N_rays_band);
    if N_rays_remaining > 0
        remaining_Q_band = Q_emit_tot_band-N_rays_band*power_per_ray; % Get remaining emissive power
        if any(remaining_Q_band(:)<0)
            remaining_Q_band(remaining_Q_band<0) = 0;
        end
        CDF_emit_band_rem = [0;cumsum(remaining_Q_band/(sum(remaining_Q_band)))]; % CDF of remaining emissions from voxels in band
        N_rays_band = N_rays_band + histcounts(rand(N_rays_remaining,1),CDF_emit_band_rem)';
    end
    
    N_rays_vx_band = round(N_rays_band.*Q_emit_vx_tot_bands./(Q_emit_tot_band)); % Not probabilistic, number of rays from voxels in each band
    N_rays_flux_band = N_rays_band-N_rays_vx_band; % Number of rays from flux in each band
    
    
    %% Initialize results arrays
    emission_counts = zeros(N_vx_tot,1); % Combined emission counts vector
    absorption_pos = cell(N_bands,1); % Absorptions positions
    self_absorption_pos = cell(N_bands,1);

    %% Monitoring
    ray_tracing_time = zeros(N_bands,1);
    ray_generation_time = zeros(N_bands,1);

    for n = 1:N_bands % Don't run this in parallel, since some bands might have very little radiation thus those workers just wait around (or may have fewer bands than workers)
        if N_rays_band(n) == 0
            continue
        end
        % Unpack parts of voxel space which are used for ray tracing for current band 
        VS_opaq = voxel_spaces{n}.opaque_voxels;
        VS_opaq_eps = voxel_spaces{n}.opaque_emissivities;
        VS_surf_norms = voxel_spaces{n}.surface_normals;
        %surf_norm_inds = voxel_spaces{n}.surface_normal_inds;
        %surf_norms_lin = voxel_spaces{n}.surface_normals_lin;
        VS_PM_kappa = voxel_spaces{n}.PM_absorption_coeffs;
        VS_nn = voxel_spaces{n}.refractive_indexes;
        fluxes = voxel_spaces{n}.fluxes;

        %% Handle all Emissions
        tic
        if N_rays_vx_band(n) > 0
            inds_vx_band = inds_vx_bands{n};  
            Q_emit_filter = VS_Q_emit_bands{n}(inds_vx_band); % Just take emitting voxels (becomes a 1D array)
            emission_counts_band = floor(N_rays_vx_band(n)*Q_emit_filter/Q_emit_vx_tot_bands(n)); % Assign rays deterministically
            N_rays_remaining_band = N_rays_vx_band(n)-sum(emission_counts_band);
            if N_rays_remaining_band > 0
                remaining_Q_emit = Q_emit_filter-emission_counts_band*power_per_ray; % Get remaining emissive power
                if any(remaining_Q_emit(:)<0)
                    remaining_Q_emit(remaining_Q_emit<0) = 0;
                end
                CDF_emit_vx = [0;cumsum(remaining_Q_emit/(sum(remaining_Q_emit)))]; % CDF of remaining emissions from voxels in band
                CDF_emit_vx = CDF_emit_vx/CDF_emit_vx(end);
                rand_vals = rand(N_rays_remaining_band,1);
                emission_counts_band = emission_counts_band + histcounts(rand_vals,CDF_emit_vx)'; % Assign remaining emissions probablistically
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
        absorption_flux_rays = []; % If flux ray is generated within an opaque voxel, it should be instantly absorbed without ray tracing
        if N_rays_flux_band(n) > 0
            N_flux = length(fluxes);
            N_rays_flux = zeros(N_flux,1); 
            rays_flux_pos = [];
            rays_flux_dir = [];
            for ii = 1:N_flux
                if ii < N_flux
                    N_rays_flux(ii) = round(N_rays_flux_band(n)*fluxes{ii}.power/Q_flux_band(n));
                else
                    N_rays_flux(ii) = N_rays_flux_band(n)-sum(N_rays_flux);
                end
                [rays_flux_pos_cur,rays_flux_dir_cur] = fluxes{ii}.GenerateRays(N_rays_flux(ii));
                inds_pos_flux = (min(floor(rays_flux_pos_cur)+1,size_VS));
                rays_pos_lin = inds_pos_flux(:,1) + (inds_pos_flux(:,2)-1)*size_VS(1) + (inds_pos_flux(:,3)-1)*size_VS(1)*size_VS(2);
                absorbed_rays_rows = VS_opaq(rays_pos_lin); % True for rays which are generated within an opaque voxel
                absorption_flux_rays = [absorption_flux_rays;inds_pos_flux(absorbed_rays_rows,:)]; %#ok<AGROW>
                rays_flux_pos = [rays_flux_pos;rays_flux_pos_cur(~absorbed_rays_rows,:)]; %#ok<AGROW> % Rather than storing as cell and using cell2mat, we just grow it since we expect the number of fluxes to usually be small
                rays_flux_dir = [rays_flux_dir;rays_flux_dir_cur(~absorbed_rays_rows,:)]; %#ok<AGROW>
            end
            rays_pos = [rays_pos;rays_flux_pos]; %#ok<AGROW>
            rays_dir = [rays_dir;rays_flux_dir]; %#ok<AGROW>
        end
        %if any(size_VS==1) % This makes traversal faster should we have a 2D or 1D domain since we don't have to constantly bounce off the boundary
        %    rays_dir(:,size_VS==1) = 0;
        %    rays_dir = rays_dir./vecnorm(rays_dir,2,2);
        %end

        ray_generation_time(n) = toc;
        tic
        [absorption_pos_cur,events] = traverseRays(rays_pos,rays_dir,VS_opaq,VS_opaq_eps,VS_surf_norms,VS_PM_kappa,VS_nn,reflective_BCs,Vxyz,size_VS);
        ray_tracing_time(n) = toc;
        try
        self_absorption_pos{n} = absorption_pos_cur(events(1:N_rays_vx_band(n))==4,:); % Self absorptions only considered if the ray is emitted from a voxel, not if emitted from an Flux object
        catch e
        debug = 0
        end
        absorption_pos{n} = [absorption_pos_cur((events==2) | (events ==3),:);absorption_pos_cur((N_rays_vx_band(n)+1:end)==4,:);absorption_flux_rays];
        
    end % Iterating over spectral bands
    absorption_pos = cell2mat(absorption_pos);
    if isempty(absorption_pos>0)
        absorption_counts = zeros(N_vx_tot,1);
    else
        absorption_counts = histcounts(absorption_pos(:,1) + ...
                                   (absorption_pos(:,2) - 1)*size_VS(1) + ...
                                   (absorption_pos(:,3) - 1)*size_VS(1)*size_VS(2), ... # This is faster version of sub2ind
                                    1:(N_vx_tot+1))';
    end
    self_absorption_pos = cell2mat(self_absorption_pos);
    if isempty(self_absorption_pos)
        self_absorption_counts = zeros(N_vx_tot,1);
    else
        self_absorption_counts = histcounts(self_absorption_pos(:,1) + ...
                                   (self_absorption_pos(:,2) - 1)*size_VS(1) + ...
                                   (self_absorption_pos(:,3) - 1)*size_VS(1)*size_VS(2), ... # This is faster version of sub2ind
                                    1:(N_vx_tot+1))';
    end
    emission_counts = emission_counts-self_absorption_counts;
    %% Calculate Delta_Q
    % Net number of absorbed rays for each voxel * power per ray = Power out of each voxel
    VS_Q_self_absorb = self_absorption_counts*power_per_ray;
    VS_Q_self_absorb = reshape(VS_Q_self_absorb,size_VS);

    VS_dQ = (absorption_counts - emission_counts)*power_per_ray; % (1D double) [W]
    VS_dQ = reshape(VS_dQ,size_VS); % (3D double) [W]: reshape to 3D voxel space

    VS_Q_emit_no_self = VS_Q_emit - VS_Q_self_absorb; % Remove self-absorptions from emissive power -> this is relavent for coupling with conduction 
    
    total_time = toc(outer_timer);
    if any(contains({'concise','verbose'},output_mode))
        fprintf("Heat Flows: total time = %0.1f s    ray generation = %0.1f s   ray tracing = %0.1f s\n", ... 
                total_time,sum(ray_generation_time),sum(ray_tracing_time));
    end
end