clear
close all
clc
format shortG
%% Parallel Plates with Participating Medium inbetween and conduction
% Cooling of Cylindrical Glass Disk - Lee and Viskanta 1998, "Transient conductive-radiative cooling of an optical
% quality glass disk", Schott BK7 data taken from https://refractiveindex.info/database/data/glass/schott/N-BK7.yml
% Note that although the refractive indexes match the values shown in the paper, the extinction coefficients do not.
% In this version we attempt to exploit symmetry in all axes to improve resolution
if isempty(gcp('nocreate'))
    parpool("Threads",8); % Threads profile reduces memory usage and also makes the pool start instantly
end
total_tic = tic;
%% Params
% Cylinder longitudinal axis is Z-direction
% Define cylinder radius and height in voxels 
params.R_vx = 30; %
params.L_vx = 40; % 
params.N_rays = 1e6; % Number of rays per step
params.time_steps = [0.1];  % Array of time steps increasing in size 
params.time_step_changes = [100]; % Tells us when to go to the next time step or (for the last value) end the simulation
params.use_internal_itr = false; % Use internal iteration on each time step > 0.1 (see Transient Cond Rad)
params.ns = 3; % Neighbourhood size for normals (only if use_manual_surface_area is false)
params.use_manual_surface_area = false;
params.spectral_band_edges = [linspace(0.5,5,8),1e5];
params.kappa_source = 'abs_coeff_paper.txt'; % 'extinction_coeff_data.txt'
params.simulation_name = 'GlassDiskCoolingSymmetric';

%% Unpack params -- seems pointless but writing the code in this way makes it very easy to functionalize
R_vx = params.R_vx;
L_vx = params.L_vx;
N_rays = params.N_rays;
time_steps = params.time_steps;
time_step_changes = params.time_step_changes;
use_internal_itr = params.use_internal_itr;
use_manual_surface_area = params.use_manual_surface_area;
ns = params.ns;
spectral_band_edges = params.spectral_band_edges;
kappa_source = params.kappa_source;
simulation_name = params.simulation_name;

%% Problem Definition

% Physical Dimensions of Disk (pp 8)
R = 100/1000; % [m]: Radius of Disk
L = 50/1000; % [m]: Thickness of Disk

% Properties of disk
density = 2514.8; % [kg/m^3]: Disk density
thermal_conductivity = 1.672; % [W /(kg K)]: Thermal conductivity
specific_heat = 1239.6; % [J /(kg K)]: Specific heat capacity
eps = 0.9; % [-]: Emissivity for lambdas > 5 um
T_disk_initial = 600+273.15; % [K]
%T_disk_initial = 1000; 

% Other
T_wall = 20+273.15; % [K]
T_surroundings = 20+273.15; % [K]:
%T_wall = 0;
%T_surroundings = 0;

% Convection
%h_top = 11.2; % [W/(m^2 K)]:
%h_bottom = 5.6; % [W/(m^2 K)]:
%h_vert = 10.7; %[W/(m^2 K)]:

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K^4];

%% Derived Parameters
% We will pad the disk in all directions with ns voxels of vacuum, followed by 1 voxel of opaque bounding surface
N_bands = max(length(spectral_band_edges)-1,1);
size_VS = [R_vx,R_vx,L_vx]+2; % Plus 2 allows for vacuum gap and then black wall.
vx_scale = [R/R_vx,R/R_vx,L/2/L_vx]; % [m/vx], note L is halfed due to z-axis symmetry
% total time steps
total_time_steps = sum([time_step_changes(1),diff(time_step_changes)]./time_steps) % Just out of interest for planning simulation time (doesn't actually do anything)


%% Open project
cur_folder = pwd;
path_parts = strsplit(cur_folder,filesep);
head = path_parts{end};
while head ~= "FIVER"
    path_parts = path_parts(1:end-1);
    head = path_parts{end};
end
FIVER_proj_path = fullfile(strjoin(path_parts,filesep),"FIVER.prj");
if isempty(matlab.project.rootProject)
    open(FIVER_proj_path);
    cd(cur_folder); % change back to TES folder
end
% Refractive index function defined at end of script!

%% Get Absorption Coefficient Data
table_data = readtable(fullfile(cur_folder,kappa_source),'VariableNamingRule','preserve');

if strcmp(kappa_source,'extinction_coeff_data.txt')
    lambdas_table = table_data{:,1}; % [um]
    kappas_table = 4*pi*table_data{:,2}./(lambdas_table*10^(-6)); % [1/m] Convert to absorption coeffs (Modest Eq 2.43)
    additional_lambdas = [2.8;4.3;5];
    additional_kappas = [10;10;15]*100; % Approximate values from Fig. 2 of Lee paper -> convert to 1/m from 1/cm
    lambdas_table = [lambdas_table;additional_lambdas];
    kappas_table = [kappas_table;additional_kappas];
elseif strcmp(kappa_source,'abs_coeff_paper.txt')
    lambdas_table = table_data{:,1}; % [um]
    kappas_table = table_data{:,2}*100; % [1/m] Convert from 1/cm to 1/m
end

%% Generate Voxel Space
voxel_spaces = cell(N_bands,1);

VS_disk = false(R_vx,R_vx,L_vx);
center = 0; % Defining quarter slice of cylinder due to triaxial symmetry
A_norm = [vx_scale(2)*vx_scale(3),vx_scale(1)*vx_scale(3),vx_scale(1)*vx_scale(2)];
VS_norms_manual = cell(size_VS);
VS_surf_areas_manual = zeros(size_VS);
for ii = 1:R_vx
    for jj = 1:R_vx
        if sqrt((ii-0.5)^2+(jj-0.5)^2) <= R_vx % if center of voxel is within circle
            VS_disk(ii,jj,1:L_vx) = true;
            VS_norms_manual{ii,jj,L_vx} = [0,0,1];
            VS_norms_manual{ii,jj,L_vx+1} = [0,0,-1];
            VS_surf_areas_manual(ii,jj,L_vx) = A_norm(3);
            if sqrt((ii+0.5)^2+(jj-0.5)^2) > R_vx 
               cur_norm = [ii-0.5,jj-0.5,0]/norm([ii-0.5,jj-0.5,0]);
               
               VS_surf_areas_manual(ii,jj,1:L_vx) = 1/max(abs(cur_norm)./A_norm);
               VS_norms_manual{ii,jj,L_vx} = cur_norm/2 + [0,0,1]/2; % average of edge vector and top surface vector;
               for kk = 1:(L_vx-1)
                VS_norms_manual{ii,jj,kk} = cur_norm;
                VS_norms_manual{ii+1,jj,kk} = -[ii+0.5,jj-0.5,0]/norm([ii+0.5,jj-0.5,0]);
               end
            end
            if sqrt((ii-0.5)^2+(jj+0.5)^2) > R_vx % If it's on the boundary
               cur_norm = [ii-0.5,jj-0.5,0]/norm([ii-0.5,jj-0.5,0]);
               
               VS_surf_areas_manual(ii,jj,1:L_vx) = 1/max(abs(cur_norm)./A_norm);
               VS_norms_manual{ii,jj,L_vx} = cur_norm/2 + [0,0,1]/2; % edge vector and top surface vector;
               for kk = 1:(L_vx-1)
                VS_norms_manual{ii,jj,kk} = cur_norm;
                VS_norms_manual{ii,jj+1,kk} = -[ii-0.5,jj+0.5,0]/norm([ii-0.5,jj+0.5,0]);
               end
            end
        end
        
    end
end
%VS_disk2 = padarray(VS_disk,[1,1,1],false,'post');
%[norms,surface_areas] = getNormalsAndSurfaceAreas(VS_disk2,vx_scale,ns);
SA = sum(VS_surf_areas_manual(:));
SA_analytic = (2*pi*R^2+R*pi*2*L)/8;
SA_error = ((SA-SA_analytic)/SA_analytic)*100;
fprintf("Surface area error = %0.4f %%\n",SA_error)
VS_surf_areas_manual(:,:,1:L_vx-1) = VS_surf_areas_manual(:,:,1:L_vx-1)*R*pi*2*L*(1-1/L_vx)/8/sum(VS_surf_areas_manual(:,:,1:L_vx-1),'all');
SA = sum(VS_surf_areas_manual(:));
SA_analytic = (2*pi*R^2+R*pi*2*L)/8;
SA_error = ((SA-SA_analytic)/SA_analytic)*100;
fprintf("Surface area error = %0.4f %%\n",SA_error)
for ii = 1:R_vx
    for jj = 1:R_vx
        if sqrt((ii-0.5)^2+(jj-0.5)^2) <= R_vx && (sqrt((ii+0.5)^2+(jj-0.5)^2) > R_vx || sqrt((ii-0.5)^2+(jj+0.5)^2) > R_vx) % If it's on the boundary % if center of voxel is within circle      
           VS_surf_areas_manual(ii,jj,L_vx) = VS_surf_areas_manual(ii,jj,L_vx-1);
        end
    end
end
SA = sum(VS_surf_areas_manual(:));
SA_analytic = (2*pi*R^2+R*pi*2*L)/8;
SA_error = ((SA-SA_analytic)/SA_analytic)*100;
fprintf("Surface area error = %0.4f %%\n",SA_error)
SA_missing = SA_analytic-SA;
SA_edge = sum(VS_surf_areas_manual(:,:,1),'all');
SA_top_edge_adjust = (SA_missing+SA_edge)/SA_edge;
for ii = 1:R_vx
    for jj = 1:R_vx
        if sqrt((ii-0.5)^2+(jj-0.5)^2) <= R_vx && (sqrt((ii+0.5)^2+(jj-0.5)^2) > R_vx || sqrt((ii-0.5)^2+(jj+0.5)^2) > R_vx) % If it's on the boundary % if center of voxel is within circle      
           VS_surf_areas_manual(ii,jj,L_vx) = VS_surf_areas_manual(ii,jj,L_vx)*SA_top_edge_adjust;
        end
    end
end
SA = sum(VS_surf_areas_manual(:));
SA_analytic = (2*pi*R^2+R*pi*2*L)/8;
SA_error = ((SA-SA_analytic)/SA_analytic)*100;
fprintf("Surface area error = %0.4f %%\n",SA_error)

VS_disk = padarray(VS_disk,[2,2,2],false,'post');
VS_wall = false(size_VS-1);
VS_wall = padarray(VS_wall,[1 1 1],true,'post');
[iis,jjs,zzs] = ind2sub(size_VS,find(VS_wall));
for nn1 = 1:length(iis)
    if iis(nn1) == size_VS(1) && jjs(nn1) ~= size_VS(2) && zzs(nn1) ~= size_VS(3)
        VS_norms_manual{iis(nn1),jjs(nn1),zzs(nn1)} = [-1,0,0];
        VS_surf_areas_manual(iis(nn1),jjs(nn1),zzs(nn1)) = A_norm(1);
    elseif iis(nn1) ~= size_VS(1) && jjs(nn1) == size_VS(2) && zzs(nn1) ~= size_VS(3)
        VS_norms_manual{iis(nn1),jjs(nn1),zzs(nn1)} = [0,-1,0];
        VS_surf_areas_manual(iis(nn1),jjs(nn1),zzs(nn1)) = A_norm(2);
    elseif iis(nn1) ~= size_VS(1) && jjs(nn1) ~= size_VS(2) && zzs(nn1) == size_VS(3)
        VS_norms_manual{iis(nn1),jjs(nn1),zzs(nn1)} = [0,0,-1];
        VS_surf_areas_manual(iis(nn1),jjs(nn1),zzs(nn1)) = A_norm(3);
    end
end

reflective_BCs = false(2,3); % Initialize reflective BCs
reflective_BCs(1,:) = true; % Lower boundary in all axes is reflective

nn_disk = zeros(N_bands,1);
kappa_disk = zeros(N_bands,1);

% These are the same for all the bands except the final one, where they are modified so that the disk is opaque.
VS_opaq = VS_wall;
VS_opaq_eps = double(VS_opaq); % Walls are black

for i = 1:N_bands
    voxel_space = VoxelSpace();        
    VS_nn = ones(size_VS);
    VS_PM_kappa = zeros(size_VS);

    if i == 1
        VS_kk = zeros(size_VS);
        VS_rho = zeros(size_VS);
        VS_cc = zeros(size_VS);
        
        VS_kk(VS_disk) = thermal_conductivity;
        VS_rho(VS_disk) = density;
        VS_cc(VS_disk) = specific_heat;

        voxel_space.thermal_conductivity = VS_kk;
        voxel_space.density = VS_rho;
        voxel_space.specific_heat = VS_cc;
    end

    if  spectral_band_edges(i) >= lambdas_table(1) && spectral_band_edges(i+1) <= lambdas_table(end) % Treat disk as PM
        % Get refractive index as an average of 10 values across the band (since we have an analytic expression for nn
        fun_nn = @(wavelength) get_refractive_index(wavelength);
        nn_disk(i) = integral(fun_nn,spectral_band_edges(i),spectral_band_edges(i+1))/(spectral_band_edges(i+1)-spectral_band_edges(i));
        
        % Get kappa disk by interpolating over 10 pt average
        fun_kappa = @(wavelength) interp1(lambdas_table,kappas_table,wavelength);
        kappa_disk(i) = integral(fun_kappa,spectral_band_edges(i),spectral_band_edges(i+1))/(spectral_band_edges(i+1)-spectral_band_edges(i));

        VS_nn(VS_disk) = nn_disk(i);
        VS_PM_kappa(VS_disk) = kappa_disk(i);

    else % lambda_avg > 5 -> Treat disk as opaque
        nn_disk(i) = 1; % It is already assigned though since we initialize VS_nn as ones
        VS_opaq(VS_disk) = true;
        VS_opaq_eps(VS_disk) = eps; % Assign emissivity
    end
    if i == 1 % PM bands (is the same in each PM band)
        if use_manual_surface_area
            VS_norms = VS_norms_manual;
            VS_surf_areas = VS_surf_areas_manual;
            VS_surf_areas(VS_disk) = 0;
        else
            [VS_norms,VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,vx_scale,ns,VS_nn); % 
        end     
    elseif i == N_bands %% opaque band
        if use_manual_surface_area
            VS_norms = VS_norms_manual;
            VS_surf_areas = VS_surf_areas_manual;
        else
        % Need to compute separately in case ns > 1, since the gap between the wall and the disk is only 1 voxel, and so the normal calculation can get
        % messed up.
        [VS_norms_disk, VS_surf_areas_disk] = getNormalsAndSurfaceAreas(VS_disk,vx_scale,ns); % Get surface norms;
        [VS_norms_wall,VS_surf_areas_wall] = getNormalsAndSurfaceAreas(VS_wall,vx_scale,ns);
        
        VS_norms = VS_norms_disk;
        VS_norms(cellfun('isempty',VS_norms)) = VS_norms_wall(cellfun('isempty',VS_norms));
        VS_surf_areas = VS_surf_areas_wall+VS_surf_areas_disk;
        end
    end
        
    voxel_space.opaque_voxels = VS_opaq;
    voxel_space.opaque_emissivities = VS_opaq_eps;
    voxel_space.surface_normals = VS_norms;
    voxel_space.surface_areas = VS_surf_areas;
    voxel_space.PM_absorption_coeffs = VS_PM_kappa;
    voxel_space.refractive_indexes = VS_nn;
    voxel_space.size = size_VS;
    voxel_space.voxel_scale = vx_scale;
    voxel_space.reflective_BCs = reflective_BCs;
    voxel_spaces{i} = voxel_space;
end
%% Define temperature field
VS_T_initial = zeros(size_VS); % Initialize temperature voxel space
VS_T_initial(VS_disk) = T_disk_initial;     
VS_T_initial(~VS_disk) = 0;

%VS_T_initial(VS_surroundings) = T_surroundings;

VS_T_fixed = false(size_VS); % Logical matrix of fixed temperatures
%VS_T_fixed(VS_wall) = true; % Fix wall temperature
VS_T_fixed(~VS_disk) = true; % Fix surrounding temperature;

VS_T_old = VS_T_initial;

%% Define indexes of areas of voxel space we want to track

r_center = 1;
Z_center = 1;
Z_top = L_vx;

T_center_center = VS_T_initial(r_center,r_center,Z_center);
T_top_center = VS_T_initial(r_center,r_center,Z_top);
cur_time = 0;
total_simulation_time = 0;
N_step_levels = length(time_step_changes);
ind = 0;
step_level =0;
fprintf('Cooling time: %0.1f s   Step computer time: %0.2f s   Total computer time: %0.2f s  \n',cur_time,0,total_simulation_time);
while step_level < N_step_levels
    step_level = step_level + 1;
    cur_time_step = time_steps(step_level);
    time_step_change = time_step_changes(step_level);
    
    if cur_time_step > 0.1 && use_internal_itr
        cur_internal_itr = true;
    else
        cur_internal_itr = false;
    end


    while time_step_change - cur_time > 10^(-6) % Accounts for floating point errors
        step_tic = tic;
        ind = ind+1;
        cur_time = cur_time + cur_time_step;
        VS_T_new = transientCondRad(N_rays,VS_T_old,VS_T_fixed,voxel_spaces,1,cur_time_step,cur_internal_itr,spectral_band_edges);
        T_top_center(ind+1) = VS_T_new(r_center,r_center,Z_top);
        T_center_center(ind+1) = VS_T_new(r_center,r_center,Z_center);
        time_values(ind+1) = cur_time;
        top_reference_table = readtable(fullfile(cur_folder,'Lee_Viskanta_fig_4_top_data.txt'));
        plot(top_reference_table{:,1},top_reference_table{:,2},'--','Color','k');
        hold on
        plot(time_values(1:ind+1),T_top_center(1:ind+1)-273.15);
        hold off;
        drawnow();
        
        
        VS_T_old = VS_T_new; % Update temperature field
        step_computation_time = toc(step_tic);
        total_simulation_time = toc(total_tic);
        fprintf('Cooling time: %0.1f s   Step computer time: %0.2f s   Total computer time: %0.2f s  \n',cur_time,step_computation_time,total_simulation_time);
    end
end
total_simulation_time = toc(total_tic);
params.total_simulation_time = total_simulation_time;
%% Plot Solutions:
color_scheme = load('PREC.mat').color_scheme;

center_reference_table = readtable(fullfile(cur_folder,'Lee_Viskanta_fig_4_center_data.txt'));
top_reference_table = readtable(fullfile(cur_folder,'Lee_Viskanta_fig_4_top_data.txt'));

f = figure;
line_width =  1;
hold on;
xlim([0,max(time_values)]);
ylim([500,max(T_center_center(:))-273.15]);

xlabel({"Time [s]"},'Interpreter','latex')
ylabel({"T [$^{\circ}$C]"},'Interpreter','latex')
hold on
line_style = '-';
p1_legend = plot(nan,nan,line_style,'Color',color_scheme(1),'LineWidth',line_width);
p2_legend = plot(nan,nan,line_style,'Color',color_scheme(2),'LineWidth',line_width);
p3_legend = plot(nan,nan,'--','Color','k','LineWidth',line_width);

plot(time_values,T_center_center-273.15,line_style,'Color',color_scheme(1),'LineWidth',line_width);
plot(time_values,T_top_center-273.15,line_style,'Color',color_scheme(2),'LineWidth',line_width);

plot(center_reference_table{:,1},center_reference_table{:,2},'--','Color','k');
plot(top_reference_table{:,1},top_reference_table{:,2},'--','Color','k');

legend_obj = legend([p1_legend,p2_legend,p3_legend], ...
    {'$r = 0,\; z = H/2$','$r = 0,\; z = H$','Reference Data\textsuperscript{1}'}, ...
    'Box','off', ...
    'Interpreter','latex');

set(gca,'TickLabelInterpreter','latex');
fontsize(f,'increase')
fontsize(f,'increase')
xlim([0,100])
ylim([505,605])

set(legend_obj,'Position',[0.60972, 0.6200, 0.2346, 0.1502])
box on;
set(gca,'XMinorTick','on','YminorTick','on');

T_center_ref = interp1(center_reference_table{:,1},center_reference_table{:,2},time_values,"linear",'extrap') + 273.15; % K
T_top_ref = interp1(top_reference_table{:,1},top_reference_table{:,2},time_values,"linear",'extrap') + 273.15; % K

T_center_diff = T_center_ref-T_center_center;
T_top_diff = T_top_ref-T_top_center;

RMSE_center = sqrt(mean(T_center_diff.^2))
RMSE_top = sqrt(mean(T_top_diff.^2))

% Save as 3 different file types and save params
if true % Change to false to avoid saving results
    filename = params.simulation_name;
    saveas(f,filename,'png'); % For easy viewing
    saveas(f,filename,'epsc'); % For LaTeX report
    saveas(f,filename,'fig'); % In case I want to edit it without running solver again
    save(strcat(filename,'Params.mat'),"params");
    
    params_summary = params;
    params_summary = rmfield(params_summary,'time_steps');
    params_summary = rmfield(params_summary,'time_step_changes');
    params_summary = rmfield(params_summary,'spectral_band_edges');
    params_summary.cooling_time = time_step_changes(end);
    params_summary.N_bands = N_bands;
    params_summary.RMSE_center = RMSE_center;
    params_summary.RMSE_top = RMSE_top;
    params_summary.note = "";
    writetable(struct2table(params_summary),strcat(filename,'ParamsSummary.txt'));
end

function n = get_refractive_index(lambda)
    % lambda in um
    % Can be a vector of lambdas

    %B1 = 1.03961212;
    %B2 = 2.31792344*10^(-4);
    %B3 = 1.01046945;
    %C1 = 6.0069867*10^(-3);
    %C2 = 2.00197144*10^(-2);
    %C3 = 1.03560653*10^2;

    B1 = 1.03961212;
    B2 = 2.31792344*10^(-1);
    B3 = 1.01046945;
    C1 = 6.0069867*10^(-3);
    C2 = 2.00197144*10^(-2);
    C3 = 1.03560653*10^2;

    n = sqrt(1 ...
        + B1*lambda.^2./(lambda.^2-C1) ...
        + B2*lambda.^2./(lambda.^2-C2) ...
        + B3*lambda.^2./(lambda.^2-C3));
    n = real(n); 
end



