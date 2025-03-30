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
    parpool("Threads",6); % Threads profile reduces memory usage and also makes the pool start instantly
end
total_tic = tic;
%% Params
% Cylinder longitudinal axis is Z-direction
% Define cylinder radius 
params.R_vx = 100; % Use multiple of 4 to ensure correct disk height (height is 1/2 radius, and then we have z-axis symmetry for 1/4)
params.N_rays = 4e7; % Number of rays per step
params.time_step_sizes = [0.1];  % Array of time steps increasing in size 
params.time_step_changes = [100]; % Tells us when to go to the next time step or (for the last value) end the simulation
params.use_internal_itr = false; % Use internal iteration on each time step > 0.1 (see Transient Cond Rad)
params.ns = 2; % Neighbourhood size for normals -> used for disk padding too
params.wavelength_edges = [linspace(0.3,5,6),1e6]; % Define lambda spectral bands (in um) (use 8 linear bands as in paper)
params.kappa_source = 'abs_coeff_paper.txt'; % 'extinction_coeff_data.txt'
params.simulation_name = sprintf('TransientDiskCoolingSymmetric_R%d_01',params.R_vx);

use_prev = false;

%% Unpack params -- seems pointless but writing the code in this way makes it very easy to functionalize
R_vx = params.R_vx;
N_rays = params.N_rays;
time_step_sizes = params.time_step_sizes;
time_step_changes = params.time_step_changes;
use_internal_itr = params.use_internal_itr;
ns = params.ns;
wavelength_edges = params.wavelength_edges;
kappa_source = params.kappa_source;
simulation_name = params.simulation_name;

%% Problem Definition

% Physical Dimensions of Disk (pp 8)
R = 100/1000; % [m]: Radius of Disk
L = 50/1000; % [m]: Thickness of Disk

% Properties of disk
rho = 2514.8; % [kg/m^3]: Disk density
kk = 1.672; % [W /(kg K)]: Thermal conductivity
cc = 1239.6; % [J /(kg K)]: Specific heat capacity
opaq_eps = 0.9; % [-]: Emissivity for lambdas > 5 um
T_disk_initial = 600+273.15; % [K]

% Other
T_wall = 20+273.15; % [K]
T_surroundings = 20+273.15; % [K]:

% Convection -> if I implement this!
% h_top = 11.2; % [W/(m^2 K)]:
% h_bottom = 5.6; % [W/(m^2 K)]:
% h_vert = 10.7; %[W/(m^2 K)]:

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K^4];

%% Derived Parameters
% We will pad the disk in all directions with ns voxels of vacuum, followed by 1 voxel of opaque bounding surface
N_bands = length(wavelength_edges)-1;
H_vx = round(L/R*R_vx/2); % half the Z dimension of disk in voxels based on physical dimensions and voxel diameter
size_VS = [R_vx+1+ns,R_vx+1+ns,H_vx+1+ns];
vx_scale = [R_vx, R_vx, H_vx]./size_VS; % [m/vx]

% total time steps
all_time_vals = 0;
for nn = 1:length(time_step_changes)
    N_steps_cur = (time_step_changes(nn)-all_time_vals(end))/time_step_sizes(nn);
    all_time_vals = [all_time_vals,all_time_vals(end)+(1:N_steps_cur)*time_step_sizes(nn)];
end
N_time_steps = length(all_time_vals); % Note: Includes 0
fprintf("N time steps: %d\n",N_time_steps)
%% Get Relevant Folders - For Saving Plots
% Get current folder
cur_file = mfilename('fullpath'); % This gets the filename
% Get plots folder
cur_folder_parts = split(cur_file,filesep); % Split at folder delimiters (windows or unix)
cur_folder = fullfile(cur_folder_parts{1:(end-1)}); % Remove the filename
for ii = length(cur_folder_parts):-1:1 % move backward thru folders until you find VoxelRayTracer folder
    if cur_folder_parts{ii} == "FIVER"
        root_folder = fullfile(cur_folder_parts{:});
        break;
    else
        cur_folder_parts(ii) = [];
    end
end

save_path = fullfile(cur_folder,'results',simulation_name); % does not have a file extension
VS_T_path = sprintf('%s_VS_T.mat',save_path);
data_table_path = sprintf('%s_data_table.txt',save_path);


% Refractive index function defined at end of script!
%% Generate Voxel Space
voxel_spaces = cell(N_bands,1);

VS_disk = generateQuarterCylinder(R_vx,H_vx);
VS_disk = padarray(VS_disk,[1+ns 1+ns 1+ns],false,'post'); % Pad with vacuum + wall

VS_wall = false(size_VS-1);
VS_wall = padarray(VS_wall,[1 1 1],true,'post');

VS_surroundings = logical((~VS_disk).*(~VS_wall));

reflective_BCs = false(2,3); % Initialize reflective BCs
reflective_BCs(1,:) = true; % Lower boundary specular in all dimensions (3 dimensional symmetry)

nn_disk = zeros(N_bands,1);
kappa_disk = zeros(N_bands,1);

% These are the same for all the bands except the final one, where they are modified so that the disk is opaque.
VS_opaq_PM = VS_wall;
VS_opaq_eps_PM = double(VS_opaq_PM); % Walls are black
VS_n = ones(size_VS);
VS_n(VS_disk) = 2; % Placeholder value for calculating surface normals at refractive interface
[VS_norms_PM, VS_surf_areas_PM] = getNormalsAndSurfaceAreas(VS_opaq_PM,vx_scale,ns,VS_n); % Get surface norms;
        

VS_opaq_no_PM = VS_wall;
VS_opaq_no_PM(VS_disk) = true;
VS_opaq_eps_no_PM = double(VS_opaq_no_PM)*opaq_eps;
[VS_norms_no_PM,VS_surf_areas_no_PM] = getNormalsAndSurfaceAreas(VS_opaq_no_PM,vx_scale,ns);
VS_n_no_PM = ones(size_VS);
VS_PM_kappa_no_PM = zeros(size_VS);

for ii = 1:N_bands
    voxel_space = VoxelSpace();
    band_l = wavelength_edges(ii);
    band_u = wavelength_edges(ii+1);
    if ii == 1 % only need to store these properties in the first band
        VS_kk = zeros(size_VS);
        VS_rho = zeros(size_VS);
        VS_cc = zeros(size_VS);
        
        VS_kk(VS_disk) = kk;
        VS_rho(VS_disk) = rho;
        VS_cc(VS_disk) = cc;
        
        voxel_space.thermal_conductivity = VS_kk;
        voxel_space.density = VS_rho;
        voxel_space.specific_heat = VS_cc;
    end



    if  band_u <= 5 && band_l >= 0.26  % Treat disk as PM
        % Get refractive index as an average of 10 values across the band (since we have an analytic expression for nn
        nn_disk(ii) = getMeanRefractiveIndex(band_l,band_u);%T_disk_initial);
        kappa_disk(ii) = getMeanBK7coeffViskanta(band_l,band_u);%,T_disk_initial,nn_disk(ii)); % 1/m;
        
        VS_n = ones(size_VS);
        VS_n(VS_disk) = nn_disk(ii);
        VS_PM_kappa = zeros(size_VS);
        VS_PM_kappa(VS_disk) = kappa_disk(ii); % 1/vx;

        voxel_space.surface_normals = VS_norms_PM;
        voxel_space.surface_areas = VS_surf_areas_PM;
        voxel_space.refractive_indexes = VS_n;
        voxel_space.opaque_voxels = VS_opaq_PM;
        voxel_space.opaque_emissivities = VS_opaq_eps_PM;
        voxel_space.PM_absorption_coeffs = VS_PM_kappa;
        voxel_space.refractive_indexes = VS_n;

    else % Treat disk as opaque
        nn_disk(ii) = 1;
        

        VS_opaq_PM(VS_disk) = 1;
        VS_opaq_eps_PM(VS_disk) = opaq_eps; % Assign emissivity
        voxel_space.surface_normals = VS_norms_no_PM;
        voxel_space.surface_areas = VS_surf_areas_no_PM;
        voxel_space.refractive_indexes = VS_n_no_PM;
        voxel_space.opaque_voxels = VS_opaq_no_PM;
        voxel_space.opaque_emissivities = VS_opaq_eps_no_PM;
        voxel_space.PM_absorption_coeffs = VS_PM_kappa_no_PM;
    end       
    voxel_space.size = size_VS;
    voxel_space.voxel_scale = vx_scale;
    voxel_space.reflective_BCs = reflective_BCs;
    voxel_spaces{ii} = voxel_space;
end

%% Define temperature field
VS_T_initial = zeros(size_VS); % Initialize temperature voxel space
VS_T_initial(VS_disk) = T_disk_initial;     
VS_T_initial(VS_wall) = T_wall;

VS_T_initial(VS_surroundings) = T_surroundings;

VS_T_fixed = false(size_VS); % Logical matrix of fixed temperatures
VS_T_fixed(VS_wall) = true; % Fix wall temperature
VS_T_fixed(VS_surroundings) = true; % Fix surrounding temperature;


if use_prev
    if isfile(VS_T_path)
        VS_T_prev = load(VS_T_path).VS_T_new;
        data_table = readtable(data_table_path);
    else % Doesn't exist
        use_prev = false;
    end
end
table_var_names = {'time_vals','T_top_center','T_center_center','step_comp_time','step_N_rays_e6'};
Z_top = size_VS(3)-ns-1;
if use_prev
    VS_T_initial = VS_T_prev;
    prev_time_vals = data_table.time_vals';
    all_time_vals(all_time_vals<=prev_time_vals(end)) = [];
    all_time_vals = [prev_time_vals,all_time_vals];
    N_time_steps = length(all_time_vals);
    N_vals_prev = length(prev_time_vals);
    last_time_val = prev_time_vals(end);
    T_top_center = [data_table.T_top_center;zeros(N_time_steps-N_vals_prev,1)];
    T_center_center = [data_table.T_center_center;zeros(N_time_steps-N_vals_prev,1)];
    prev_total_time = sum(data_table.step_comp_time);
    
else
    prev_time_vals = 0;
    T_top_center = zeros(N_time_steps,1);
    T_center_center = zeros(N_time_steps,1);
    T_top_center(1,:) = T_disk_initial;
    T_center_center(1,:) = T_disk_initial;
    data_table = table(0,T_disk_initial,T_disk_initial,0,0,'VariableNames',table_var_names);
    prev_total_time = 0;
end
VS_T_prev = VS_T_initial;

%% Define indexes of areas of voxel space we want to track

cur_time = 0;
total_simulation_time = 0;
N_step_levels = length(time_step_changes);
fprintf('Cooling time: %0.1f s   Step computer time: %0.2f s   Total computer time: %0.2f s  \n',cur_time,0,total_simulation_time);
prev_time = 0;


%% Prepare figure
f = figure;
ax = gca;
line_width =  1;
hold on;
xlim([0,all_time_vals(end)]);
ylim([500,max(T_center_center(:))-273.15]);
xlabel({"Time [s]"},'Interpreter','latex')
ylabel({"T [$^{\circ}$C]"},'Interpreter','latex')

color_scheme = load(fullfile(root_folder,'src','data','ColorSchemes','PREC.mat')).color_scheme;
center_reference_table = readtable(fullfile(cur_folder,'Lee_Viskanta_fig_4_center_data.txt'));
top_reference_table = readtable(fullfile(cur_folder,'Lee_Viskanta_fig_4_top_data.txt'));


line_style = '-';
p1_legend = plot(nan,nan,line_style,'Color',color_scheme(1),'LineWidth',line_width);
p2_legend = plot(nan,nan,line_style,'Color',color_scheme(2),'LineWidth',line_width);
p3_legend = plot(nan,nan,'--','Color','k','LineWidth',line_width);

p_ref1 = plot(center_reference_table{:,1},center_reference_table{:,2},'--','Color','k');
p_ref2 = plot(top_reference_table{:,1},top_reference_table{:,2},'--','Color','k');

p_cc = plot(prev_time_vals,T_center_center(1:length(prev_time_vals))-273.15,line_style,'Color',color_scheme(1),'LineWidth',line_width);
p_tc = plot(prev_time_vals,T_top_center(1:length(prev_time_vals))-273.15,line_style,'Color',color_scheme(2),'LineWidth',line_width);

legend_obj = legend([p1_legend,p2_legend,p3_legend], ...
    {'$r = 0,\; z = H/2$','$r = 0,\; z = H$','Reference Data\textsuperscript{1}'}, ...
    'Box','off', ...
    'Interpreter','latex');
set(legend_obj,'Position',[0.60972, 0.6200, 0.2346, 0.1502]);

fontsize(f,'increase')
fontsize(f,'increase')

drawnow();


%% Run simulation
for ind = 1:N_time_steps %
    
    cur_time = all_time_vals(ind);
    cur_time_step = cur_time-prev_time;
    prev_time = cur_time;
    if cur_time <= prev_time_vals(end)
        continue
    end
    step_tic = tic;
    VS_T_new = transientCondRad(N_rays,VS_T_prev,VS_T_fixed,voxel_spaces,1,cur_time_step,use_internal_itr,wavelength_edges);
    
    T_top_center(ind) = VS_T_new(1,1,Z_top);
    T_center_center(ind) = VS_T_new(1,1,1);
    
    VS_T_prev = VS_T_new; % Update temperature field
    step_computation_time = toc(step_tic);
    total_simulation_time = toc(total_tic)+prev_total_time;
    fprintf('Cooling time: %0.1f s   Step computer time: %0.2f s   Total computer time: %0.2f s  \n',cur_time,step_computation_time,total_simulation_time); 
    % Update plots
    set(ax.Children(2),'XData',all_time_vals(1:ind),'YData',T_center_center(1:ind)-273.15); % p_cc
    set(ax.Children(1),'XData',all_time_vals(1:ind),'YData',T_top_center(1:ind)-273.15); % Note new children are added to axes at index 1, so this is p_tc
    drawnow();

    data_table(end+1,:) = table(all_time_vals(ind),T_top_center(ind),T_center_center(ind),step_computation_time,round(N_rays/1e6)); % Add row to table
    
    % Save as we go!
    writetable(data_table,data_table_path);
    save(VS_T_path,'VS_T_new');
    saveas(f,save_path,'fig');
    saveas(f,save_path,'png');
    saveas(f,save_path,'epsc');
end
total_simulation_time = toc(total_tic) + prev_total_time;

function abs_coeff = getBK7coeffViskanta(wavelength,data_table)
% Input:
%   wavelength (um)
% Output:
%   abs_coef (1/m)
    if nargin == 1
        data_table = readtable('BK7_data_Viskanta.txt');
    end
    if wavelength > 5
        error('Data not available above 5 um')
    elseif wavelength<0.26
        error('Data not available below 0.26 um')
    end
    abs_coeff = interp1(data_table.wavelength,data_table.abs_coeff,wavelength)*100; % Convert from 1/cm to 1/m
end

function mean_abs_coeff = getMeanBK7coeffViskanta(wavelength1,wavelength2)
    data_table = readtable('BK7_data_Viskanta.txt');
    mean_abs_coeff = integral(@(x) getBK7coeffViskanta(x,data_table),wavelength1,wavelength2)/(wavelength2-wavelength1);

end

function f = brf(wavelength,T,n_refrac)
    wavelength_m = wavelength*1e-6; % convert from um to m
    h = 6.62607004*10^(-34);     % Planck's constant [J s]
    c = 299792458;              % speed of light in a vacuum [m/s]
    Kb = 1.38064852*10^(-23);    % Boltzmann's constant [J/K]
    f = (2*pi*h*c^2)./(n_refrac.^2*wavelength_m.^5.*(exp((h*c)./(n_refrac.*wavelength_m.*Kb.*T))-1));
end

function mean_abs_coeff = getBRF_MeanBK7coeffViskanta(wavelength1,wavelength2,T,n_refrac)
    data_table = readtable('BK7_data_Viskanta.txt');
    mean_abs_coeff = integral(@(x) getBK7coeffViskanta(x,data_table).*brf(x,T,n_refrac),wavelength1,wavelength2)/(wavelength2-wavelength1)/spectralBandPower(lamda1,lambda2,T,n_refrac);

end

function n = getRefractiveIndex(lambda)
    % lambda in um
    % Can be a vector of lambdas

    B1 = 1.03961212;
    B2 = 2.31792344*10^(-4);
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

function mean_n = getMeanRefractiveIndex(lambda1,lambda2)
    mean_n = integral(@(x) getRefractiveIndex(x),lambda1,lambda2)/(lambda2-lambda1);
end

function mean_n = getBRF_MeanRefractiveIndex(lambda1,lambda2,T)
    n_refrac = 1;
    for ii = 1:3 % Iterate a few times
        f = @(x) getRefractiveIndex(x).*brf(x,T,n_refrac);
        n_refrac = integral(f,lambda1,lambda2)/spectralBandPower(lambda1,lambda2,T,n_refrac); % convert wavelength to m in brf
    end
    mean_n = n_refrac;
end

function VS = generateQuarterCylinder(radius,height)
    VS = false(radius,radius,height);
    r_sq = radius^2;
    for i = 1:radius
        x_sq = (i-0.5)^2;
        for j = 1:radius
            y_sq = (j-0.5)^2;
            if x_sq + y_sq <= r_sq
                VS(i,j,:) = true;
            end
        end
    end
end
