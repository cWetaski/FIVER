% Author: Charles Wetaski

%% Coupled Parallel Black Plates with Spectrally Dependent Participating Medium at Equilibrium
% parallel, infinite plates (note that this is really a 1D problem)
% This follows:
% "Interaction of Heat Transfer by Conduction and Radiation in a Nongray Planar Medium" by Crosbie et.al, 1971


clear all
close all
clc
if isempty(gcp('nocreate'))
    parpool('Threads',4);
end
outer_tic = tic;
%% Params
% Plates in YZ direction, separated by distance X
% Voxel Space Dimensions
X = 100;
Y = 2; % Y boundary specularly reflective
Z = 2; % Z boundary specularly reflective

% Spectral bands
nu_bar_bands = {3,3,3,3}; % Nondimensional frequency nu_bar = h*nu/(kb*T2); where kb is boltzmann constant 
spectral_absorption_coeffs = {[1,0],[0,1],[1,0],[0,1]}; % scaling factor of optical depth in each band



% Nondimensional conduction-to-radiation interaction parameter == k*beta/(4*sigma*T_2^4) (Different from modest)
N_param = [0,0,0.05,0.05];
N_cases = length(N_param);
N_param_strings = string(N_param);

% Free parameters
T2 = 500; % [K]: Temperature of plate 2 (arbitrary)

% Initial nondimension temperature in medium
theta_0 = 0.75;

% Number of rays per iteration
N_rays_set = [2e5,5e5,1e6,2e6];
pure_rad_scale = 4; % Multiply last number of rays by this number in pure radiation case since it's much noisier;

% Convergence
max_itr = 50;
N_prev_itr = 5;

N_Psi = 1;
N_rays_Psi = 1e6;

vx_scale = [1,100,100];

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K^4]; Stefan Boltzmann
planck = 6.62607004*10^(-34);     % Planck's constant [J s]
c = 299792458;              % speed of light in a vacuum [m/s]
kb = 1.38064852*10^(-23);    % Boltzmann's constant [J/K]

tau_L = 1;

T1 = T2/2; % [K]: Results Crosbie are for this case
eps1 = 1; % emissivity of plate 1; results from Modest are just for black plates.
eps2 = 1; % emissivity of plate 2;

VS_T_equil = cell(N_cases,1);
nondim_T = cell(N_cases,1); % Theta, in crosbie paper
Psi = zeros(N_cases,1); % Nondimensional heat flux
Psi_std = zeros(N_cases,1); 

for i = 1:N_cases
    %% Derived Parameters
    size_VS = [X+2,Y,Z]; %  +2 because each plate is 1 vx thick
    N_bands = length(nu_bar_bands{i})+1;
    if isempty(nu_bar_bands{i})
        spectral_bands_edges = [];
    else
        spectral_bands_edges = [0,flip(planck*c./(nu_bar_bands{i}*kb*T2))*10^(6),1e6]; % [um]: flip order since increasing nu -> decreasing lambda 
    end
    cur_spectral_absorption_coeffs = flip(spectral_absorption_coeffs{i}); % flip order to correspond with lambda bands correctly
    PM_kappa = cur_spectral_absorption_coeffs*tau_L/(X*vx_scale(1)); % [1/m]:
    
    beta = tau_L/(X*vx_scale(1)); % 1/m
    thermal_conductivity = N_param(i)*4*sigma*T2^3/beta; % [W/(m-K)]: By definition of N
    cp = 1;
    rho = 1;
    if N_param(i) == 0
        if length(N_rays_set)>3
            N_rays = N_rays_set(end-3:end);
        else
            N_rays = N_rays_set;
        end
        N_rays(end) = N_rays(end)*pure_rad_scale; % More rays for pure radiation case
        C_relax = 1;
        C_converge = 1.5;
    elseif N_param(i) < 0.1
        C_relax = 0.7; % Underrelaxing is necessary for this case
        C_converge = 1.1;
    elseif N_param(i) > 10
        N_rays = N_rays_set(1); % Ray tracing doesn't actually matter here
        C_relax = 1;
        C_converge = 1.1;
    else
        N_rays = N_rays_set;
        C_relax = 0.7;
        C_converge = 1.1;
    end

    %% Generate Voxel Spaces
    voxel_spaces = cell(N_bands,1);
    for j = 1:N_bands
        VS_plate1 = false(size_VS);
        VS_plate1(1,:,:) = 1;
        VS_plate2 = false(size_VS);
        VS_plate2(end,:,:) = 1;
        
        VS_opaq = logical(VS_plate1 + VS_plate2); % Join the two plates into 1 opaque voxel space;
        VS_opaq_eps = double(VS_opaq); % Both plates are black bodies
        
        [VS_surf_norms, VS_surf_areas, ~] = getNormalsAndSurfaceAreas(VS_opaq,vx_scale,1); % Get surface normals and areas
        
        reflective_BCs = false(2,3); % Initialize reflective BCs
        reflective_BCs(:,2:3) = 1; % Y and Z boundaries are reflective 
        
        VS_PM_kappa = zeros(size_VS); % Initialize PM voxel space
        VS_PM_kappa(2:(end-1),:,:) = PM_kappa(j); % only thing that changes from band to band
        VS_nn = ones(size_VS); % refractive indexes can all be 1
        

        
        % Store as voxel space object
        voxel_space = VoxelSpace();
        voxel_space.opaque_voxels = VS_opaq;
        voxel_space.opaque_emissivities = VS_opaq_eps;
        voxel_space.surface_normals = VS_surf_norms;
        voxel_space.surface_areas = VS_surf_areas;
        voxel_space.PM_absorption_coeffs = VS_PM_kappa;
        voxel_space.refractive_indexes = VS_nn;
        voxel_space.thermal_conductivity = thermal_conductivity;
        voxel_space.specific_heat = 1;
        voxel_space.density = 1;
        voxel_space.size = size_VS;
        voxel_space.voxel_scale = vx_scale;
        voxel_space.reflective_BCs = reflective_BCs;
        voxel_spaces{j} = voxel_space;
    end

    %% Generate Temperature Field
    VS_T = zeros(size_VS); % Initialize temperature voxel space
    VS_T(1,:,:) = T1;
    VS_T(end,:,:) = T2;
    VS_T(2:(end-1),:,:) = theta_0*T2; % Initial temperature in medium.
    
    %% Fixed temperatures for equilibrium solver
    VS_T_fixed = false(size_VS); % Logical matrix of fixed temperatures
    VS_T_fixed(1,:,:) = 1;
    VS_T_fixed(end,:,:) = 1;
    
    n_step = 0;
    VS_T_old = VS_T;
    
    fprintf("Solving %d of %d    time = %0.3f \n",i,N_cases,toc(outer_tic))
    [VS_T_equil{i}, ~, ~, count_itr(i)] = equilibriumCondRad(N_rays,VS_T_old,VS_T_fixed,voxel_spaces, ... 
        "SpectralBandEdges",spectral_bands_edges,"OutputMode","verbose", ...
        "ConvergenceConstant",C_converge,'NPreviousIterations',N_prev_itr,'MaxIterations',max_itr, ...
        "RelaxationConstant",C_relax);
    VS_T_new = VS_T_equil{i};

    nondim_T{i} = mean(VS_T_equil{i},[2 3])/T2; % theta in Crosbie paper
    d_theta_d_tau = (nondim_T{i}(2)-nondim_T{i}(1))/(tau_L/X); % Approximation of derivative at boundary
    Psi_test = zeros(N_Psi,1);
    for j = 1:N_Psi % Calculate radiative heat flux N_psi times
        VS_dQ = radiativeHeatFlowsMC(N_rays_Psi,VS_T_equil{i},voxel_spaces, ...
            'SpectralBandEdges',spectral_bands_edges,'OutputMode','quiet');
        
        nondim_radiative_flux = -mean(VS_dQ/(vx_scale(2)*vx_scale(3)),[2,3])/(sigma*T2^4); % big F+ in Crosbie Paper, Dividing by vx_scale^2 to convert from flux divergence to flux 
        Psi_test(j) = 4*N_param(i)*d_theta_d_tau-nondim_radiative_flux(1); % Eq 33 in paper, Note that VS_dQ = radiative flux in
        
    end
    Psi(i) = mean(Psi_test);
    Psi_std(i) = std(Psi_test); 
    
    nondim_T{i} = nondim_T{i}(2:(end-1)); % Don't need endpoint values since they are fixed (T1 and T2) and make the plot weird
    %% Calculate heat flux
end
Psi(Psi==0) = [];
Psi_std(Psi_std==0) = [];
total_simulation_time = toc(outer_tic);


%% Plot Solutions - This section is a mess since I quickly wanted to change the plot and just used a bunch of random continues
color_scheme = load('PREC.mat').color_scheme;
line_width = 1;

f = figure;
firstax = axes(f);
hold on;
xlim([0,1]);
ylim([0.5,1]);
set(gca,'LooseInset',get(gca,'TightInset'));
f.Position = [f.Position(1),f.Position(2),484,396];

xlabel({"Nondimensional location","$x/L$"},'Interpreter','latex')
ylabel({"Nondimensional temperature", "$T/T_1$"},'Interpreter','latex')
yticks([0.5,0.6,0.7,0.8,0.9,1.0])
set(gca,'TickLabelInterpreter','latex');
dX = 1/X;
x_vals = linspace(dX/2,1-dX/2,X);

legend_str = strings(N_cases,1);
j = 0;

% Load ref data
ref_interp = zeros(4,X);
T_error = zeros(4,X);
ref_str = ["CrosbieModelA_N0.txt","CrosbieModelB_N0.txt","CrosbieModelA_N005.txt","CrosbieModelB_N005.txt"];
for i = 1:4
    ref_data = readtable(ref_str(i));
    ref_interp(i,:) = interp1(ref_data{:,'Var1'},ref_data{:,'Var2'},x_vals,'linear','extrap');
    T_error(i,:) = (ref_interp(i,:)'-nondim_T{i})./(nondim_T{i})*100;
    max_T_error(i) = max(abs(T_error(i,:)),[],'all');
    RMSE_T(i) = sqrt(mean(T_error(i,:).^2,'all'));
end
Psi_ref = [0.7860;0.6703;0.9051;0.8029];
Psi_errors = (Psi-Psi_ref)./Psi_ref*100;
fprintf("Computation time: %0.2f seconds\n",total_simulation_time)
fprintf("RMSE errors: %0.2f %%, %0.2f %%, %0.2f %%, %0.2f %% \n",RMSE_T(1),RMSE_T(2),RMSE_T(3),RMSE_T(4));
fprintf("max T errors: %0.2f %%, %0.2f %%, %0.2f %%, %0.2f %% \n",max_T_error(1),max_T_error(2),max_T_error(3),max_T_error(4));

data_table = table();
% Plot results
plot_count = 0;
j = 0;
for i = 1:N_cases
    if isequal(spectral_absorption_coeffs{i},[1,0])
        line_style = "-";
        plot_color = color_scheme(1);
        j = j + 1;
        legend_str(j) = "Model A";
        legendp1(j) = plot(nan,nan,'-','Color',plot_color,'LineWidth',line_width,'Parent',firstax);
    elseif isequal(spectral_absorption_coeffs{i},[0,1])
        line_style = "-";
        j = j+1;
        plot_color = color_scheme(2);
        legend_str(j) = "Model B";
        legendp1(j) = plot(nan,nan,'-','Color',plot_color,'LineWidth',line_width,'Parent',firstax);
    end
    plot_count = plot_count+1;
    plot(1-x_vals,nondim_T{i},'Color',plot_color,'LineStyle',line_style,'LineWidth',1);
    data_table(strcat('x_vals_',legend_str,'_','N=',N_param_strings(i))) = interp1;
    data_table(strcat('x_vals_',legend_str,'_','N=',N_param_strings(i))) = interp1
end
    j = j+1;
    line_style = "--";
    legend_str(j) = "Reference Data";
    legendp1(j) = plot(nan,nan,'--','Color','k','LineWidth',line_width,'Parent',firstax);
    % if N_param(i) == 0
    %     plot_color = color_scheme(1);
    % elseif N_param(i) == 0.05
    %     plot_color = color_scheme(2);
    % else
    %     continue
    %     plot_color = color_scheme(3);
    % end
    nd

% Plot Reference data
for i = 1:4
    line_style = '--';
    plot(1-x_vals,ref_interp(i,:),'Color','k','LineStyle',line_style,'LineWidth',line_width);
end

% Legend
[~,unique_inds] = unique(vertcat(legendp1.Color),'rows','stable');
legendp1 = legendp1(unique_inds);

legend1 = legend(firstax,legendp1,cellstr(unique(legend_str,'stable')), ...
    'Location','northeast', ...
    'Box','off', ...
    'Interpreter','latex');

%% Labelling the plot - This part is so messed up since I completely changed everything and mirrored the x-axis it's basically encrypted!
% Won't work properly unless following pattern of N_param and spectral_absorption_coeffs is used
% N_param = [0,0,0,0.05,0.05,0.05,100];
% spectral_absorption_coeffs = {1,[1,0],[0,1],1,[1,0],[0,1],1};

unique_N = unique(N_param,'stable');
N_strings = ["$N = 0$","$N = 0.05$","$N = 100$"];


x_begin = [0.73,0.55,0.43]; % Text x position
y_begin = [0.75,0.65,0.64]; % Text y position

x_end(1,1:3) = x_begin(1) - 0.22; % For N = 0
x_end(1,3) = x_end(1,3) - 0.05; % For model B
x_end(2,1:3) = x_begin(2) - 0.22; % For N = 0.05
x_end(2,3) = x_end(2,3) - 0.03; % For model B
x_end(3,1:3) = x_begin(3)-0.07; % For N = 1

x_begin_adj = [-0.115,-0.16,-0.01]; % Adjustment for line origins
y_begin_adj = [0.015,0.015,0.005]; 
for i = 1:3
    for count = 1:3
        [~,x_end_index(i,count)] = min(abs(x_vals-x_end(i,count))); % Get corresponding x_index
    end
end

for i = 1:length(unique_N)
    count = 0;
    label_bool = true;
    for j = 1:length(N_param)
        if N_param(j) ~= unique_N(i)
            continue;
        end
        count = count + 1;
        y_end = nondim_T{j}(x_end_index(i,count)); % Get y-position of line
        if label_bool
            N_label(i) = text(1-x_begin(i),y_begin(i),N_strings(i),'Interpreter','latex'); % Label
            label_bool = false;
        end
        hline(j) = annotation('line'); % Line callouts
        hline(j).Parent = f.CurrentAxes; % Position in units of data
        hline(j).X = [1-(x_begin(i)+x_begin_adj(i)), 1-x_end(i,count)];
        hline(j).Y = [y_begin(i)+y_begin_adj(i), y_end];
    end
end
box on
set(gca,'XMinorTick','on','YMinorTick','on')
fontsize(f,'increase')
fontsize(f,'increase')