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
R_vx = 30; %
L_vx = 40; % 
N_rays = 5e5; % Number of rays per step
time_step_sizes = [0.05];  % [seconds] Array of time steps increasing in size 
time_step_changes = [100]; % [seconds] Tells us when to go to the next time step or (for the last value) end the simulation
use_internal_itr = false; % Use internal iteration on each time step > 0.1 (see Transient Cond Rad)
ns = 1; % Neighbourhood size for normals
spectral_band_edges = [linspace(0.5,5,8),1e5]; % [um] wavelengths of spectral bands
kappa_source = 'BK7_absorptivity.txt'; % Path to file with absorptivity data
simulation_name = 'TransientDiskCooling';

%% Problem Definition
% Physical Dimensions of Disk
R = 100/1000; % [m]: Radius of Disk
L = 50/1000; % [m]: Thickness of Disk

% Properties of disk
density = 2514.8; % [kg/m^3]: Disk density
thermal_conductivity = 1.672; % [W /(kg K)]: Thermal conductivity
specific_heat = 1239.6; % [J /(kg K)]: Specific heat capacity
eps = 0.9; % [-]: Emissivity for lambdas > 5 um
T_disk_initial = 600+273.15; % [K]
% Other
T_wall = 20+273.15; % [K]

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K^4];

%% Derived Parameters
% We will pad the disk in all directions with ns voxels of vacuum, followed by 1 voxel of opaque bounding surface
N_bands = max(length(spectral_band_edges)-1,1);
size_VS = [R_vx,R_vx,L_vx]+2; % Plus 2 allows for vacuum gap and then black wall.
vx_scale = [R/R_vx,R/R_vx,L/2/L_vx]; % [m/vx], note L is halved due to z-axis symmetry
total_time_steps = sum([time_step_changes(1),diff(time_step_changes)]./time_step_sizes); 


%% Open project file if it's not already open
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
lambdas_table = table_data{:,1}; % [um]
absorptivity_table = table_data{:,2}*100; % [1/m] Convert from 1/cm to 1/m

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
        end 
    end
end
VS_disk = padarray(VS_disk,[2,2,2],false,'post'); % Pad with empty voxels
VS_wall = false(size_VS-1);
VS_wall = padarray(VS_wall,[1 1 1],true,'post'); % Voxel space of enclosure black walls

reflective_BCs = false(2,3); % Initialize reflective BCs
reflective_BCs(1,:) = true; % Lower boundary in all axes is reflective
nn_disk = zeros(N_bands,1);
kappa_disk = zeros(N_bands,1);
% These are the same for all the bands except the final one, where they are modified so that the disk is opaque.
VS_opaq = VS_wall;
VS_opaq_eps = double(VS_opaq); % Walls are black
for i = 1:N_bands % Iterate over spectral bands to create spectral voxel space
    voxel_space = VoxelSpace();        
    VS_nn = ones(size_VS);
    VS_PM_kappa = zeros(size_VS);
    if i == 1
        VS_thermal_conductivity = zeros(size_VS);
        VS_density = zeros(size_VS);
        VS_specific_heat = zeros(size_VS);
        VS_thermal_conductivity(VS_disk) = thermal_conductivity;
        VS_density(VS_disk) = density;
        VS_specific_heat(VS_disk) = specific_heat;
        voxel_space.thermal_conductivity = VS_thermal_conductivity;
        voxel_space.density = VS_density;
        voxel_space.specific_heat = VS_specific_heat;
    end
    if  spectral_band_edges(i) >= lambdas_table(1) && spectral_band_edges(i+1) <= lambdas_table(end) % Treat disk as PM
        % Get refractive index as an average of 10 values across the band (since we have an analytic expression for nn
        fun_nn = @(wavelength) get_refractive_index(wavelength);
        nn_disk(i) = integral(fun_nn,spectral_band_edges(i),spectral_band_edges(i+1))/(spectral_band_edges(i+1)-spectral_band_edges(i));
        
        % Get kappa disk by interpolating over 10 pt average
        fun_kappa = @(wavelength) interp1(lambdas_table,absorptivity_table,wavelength);
        kappa_disk(i) = integral(fun_kappa,spectral_band_edges(i),spectral_band_edges(i+1))/(spectral_band_edges(i+1)-spectral_band_edges(i));

        VS_nn(VS_disk) = nn_disk(i);
        VS_PM_kappa(VS_disk) = kappa_disk(i);

    else % lambda_avg > 5 -> Treat disk as opaque
        nn_disk(i) = 1; % It is already assigned though since we initialize VS_nn as ones
        VS_opaq(VS_disk) = true;
        VS_opaq_eps(VS_disk) = eps; % Assign emissivity
    end
    if i == 1 % PM bands (is the same in each PM band)
        [VS_norms,VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,vx_scale,ns,VS_nn); % 
    elseif i == N_bands %% opaque band
        % Need to compute separately in case ns > 1, since the gap between the wall and the disk is only 1 voxel, and so the normal calculation can get
        % messed up.
        [VS_norms_disk, VS_surf_areas_disk] = getNormalsAndSurfaceAreas(VS_disk,vx_scale,ns); % Get surface norms;
        [VS_norms_wall,VS_surf_areas_wall] = getNormalsAndSurfaceAreas(VS_wall,vx_scale,ns);
        
        VS_norms = VS_norms_disk;
        VS_norms(cellfun('isempty',VS_norms)) = VS_norms_wall(cellfun('isempty',VS_norms)); % Combine wall surface norms and disk surface norms into one cell array
        VS_surf_areas = VS_surf_areas_wall+VS_surf_areas_disk; % Combine areas
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

VS_T_fixed = false(size_VS); % Logical matrix of fixed temperatures
VS_T_fixed(~VS_disk) = true; % Fix surrounding temperature;

VS_T_old = VS_T_initial;

%% Define indexes of locations of voxel space we want to track transient temperature for
r_center = 1;
Z_center = 1;
Z_top = L_vx;

T_center_center = zeros(total_time_steps+1,1);
T_top_center = zeros(total_time_steps+1,1);
time_values = zeros(total_time_steps+1,1);

T_center_center(1) = VS_T_initial(r_center,r_center,Z_center);
T_top_center(1) = VS_T_initial(r_center,r_center,Z_top);

%% Initialize simulation
cur_time = 0;
total_simulation_time = 0;
N_step_levels = length(time_step_changes);
ind = 0;
step_level =0;

fprintf('Cooling time: %0.1f s   Step computer time: %0.2f s   Total computer time: %0.2f s  \n',cur_time,0,total_simulation_time);
figstruct = initializeFigure(VS_T_initial,VS_disk);

myVideo = VideoWriter('GlasDiskCooling','MPEG-4'); %open video file
myVideo.FrameRate = 60;  %can adjust this
open(myVideo);
frame = getframe(gcf);
writeVideo(myVideo,frame);

%% Time step loop
while step_level < N_step_levels
    step_level = step_level + 1;
    cur_time_step = time_step_sizes(step_level);
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

        %% Time step
        VS_T_new = transientCondRad(N_rays,VS_T_old,VS_T_fixed,voxel_spaces,1,cur_time_step,cur_internal_itr,spectral_band_edges);
        
        %% Update results
        T_top_center(ind+1) = VS_T_new(r_center,r_center,Z_top);
        T_center_center(ind+1) = VS_T_new(r_center,r_center,Z_center);
        time_values(ind+1) = cur_time;
        figstruct = updateFigure(figstruct,VS_T_new,time_values(1:ind+1),T_center_center(1:ind+1),T_top_center(1:ind+1));
        frame = getframe(gcf);
        writeVideo(myVideo,frame);

        VS_T_old = VS_T_new; % Update temperature field
        step_computation_time = toc(step_tic);
        total_simulation_time = toc(total_tic);
        fprintf('Cooling time: %0.1f s   Step computer time: %0.2f s   Total computer time: %0.2f s  \n',cur_time,step_computation_time,total_simulation_time);
    end
end
close(myVideo);
total_simulation_time = toc(total_tic);
params.total_simulation_time = total_simulation_time;

%% Plot Solutions:
color_scheme = load('PREC.mat').color_scheme;

% Load reference solutions
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
if false
    filename = "GlassDiskCoolingAutoSA"
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

function figstruct = initializeFigure(VS_T,VS_disk)

    center_reference_table = readtable('Lee_Viskanta_fig_4_center_data.txt');
    top_reference_table = readtable('Lee_Viskanta_fig_4_top_data.txt');
    size_VS = size(VS_T);
    f = figure;
    scattercolor = "#43A2CA";
   
    ht = tiledlayout(2,2,"TileSpacing",'compact','Padding','tight');
    f.Position = [f.Position(1)-200,f.Position(2),f.Position(3)*1.2,f.Position(4)*1.3];
    nexttile(1);
    him1 = imagesc(VS_T(:,:,1)-273.15);
    set(gca,'YDir','normal')
    set(him1,'XData',[0.5,size_VS(1)-0.5]*100/(size_VS(1)-2),'YData',[0.5,size_VS(2)-0.5]*100/(size_VS(2)-2));
    xlabel("x (mm)",'Interpreter','latex');
    ylabel("y (mm)",'Interpreter','latex')
    xlim([0,100])
    ylim([0,100])
    colorbar(); 
    set(him1,"AlphaData",VS_disk(:,:,1))
    clim([400,600]);

    hold on;
    s = scatter(0.5*100/size_VS(1),0.5*100/size_VS(2),'*','MarkerEdgeColor',scattercolor,'MarkerFaceColor',scattercolor);
    legend(s,'Sample point','Location','northeast')
    title("Center plane",'Interpreter','tex')
    

    nexttile(2);
    href2 = plot(center_reference_table{:,1},center_reference_table{:,2},'--','Color','k');
    hold on;
    hp2 = plot(nan,nan,'LineWidth',1);
    title("Sample point temperature evolution",'Interpreter','tex');
    xlim([0,100]);
    ylim([500,610]);
    xlabel({"Time [s]"},'Interpreter','latex')
    ylabel({"T [$^{\circ}$C]"},'Interpreter','latex')
    legend([hp2,href2],["FIVER","Reference data"],'Location','southwest','BackgroundAlpha',0);
    

    nexttile(3);
    set(gca,'YDir','normal')
    him3 = imagesc(VS_T(:,:,end-2)-273.15);
    set(gca,'YDir','normal')
    set(him3,'XData',[0.5,size_VS(1)-0.5]*100/(size_VS(1)-2),'YData',[0.5,size_VS(2)-0.5]*100/(size_VS(2)-2));
    xlabel("x (mm)",'Interpreter','latex')
    ylabel("y (mm)",'Interpreter','latex')
    xlim([0,100])
    ylim([0,100])
    colorbar(); 
    set(him3,"AlphaData",VS_disk(:,:,end-2))
    clim([400,600]);
    title("Top plane",'Interpreter','tex');
    hold on
    s = scatter(0.5*100/size_VS(1),0.5*100/size_VS(2),'*','MarkerEdgeColor',scattercolor,'MarkerFaceColor',scattercolor);
    legend(s,'Sample point','Location','northeast')
    

    nexttile(4);
    href4 = plot(top_reference_table{:,1},top_reference_table{:,2},'--','Color','k');
    hold on
    hp4 = plot(nan,nan,'LineWidth',1);
    title("Sample point temperature evolution",'Interpreter','tex');
    xlim([0,100]);
    ylim([500,610]);
    xlabel({"Time [s]"},'Interpreter','latex')
    ylabel({"T [$^{\circ}$C]"},'Interpreter','latex')
    legend([hp4,href4],["FIVER","Reference data"],'Location','southwest');
    
    
    drawnow();
   figstruct.f = f;
   figstruct.him1 = him1;
   figstruct.hp2 = hp2;
   figstruct.him3 = him3;
   figstruct.hp4 = hp4;
end

function figstruct = updateFigure(figstruct,VS_T,time_values,T_center_center,T_top_center)
    set(0,'CurrentFigure',figstruct.f);
    set(figstruct.him1,'CData',VS_T(:,:,1)-273.15)
    set(figstruct.hp2,'XData',time_values,'YData',T_center_center-273.15);
    set(figstruct.him3,'CData',VS_T(:,:,end-2)-273.15);
    set(figstruct.hp4,'XData',time_values,'YData',T_top_center-273.15);
    drawnow();
end



