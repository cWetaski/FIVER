clear
close all
clc
%% Parallel Plates with Participating Medium inbetween and conduction
% parallel, infinite plates (note that this is really a 1D problem)
% This reproduces the solutions from Modest example 22.1. Figure 22-1. Note that the pure radiation case (N = 0) was not
% computed since the pure radiation case is calculated for MC_Equilibrium testing, and the coupled conduction solver
% fails if you have no conduction at all. 
if isempty(gcp('nocreate'))
    parpool('Threads',6);
end
outer_tic = tic;
%% Params
% Plates in YZ direction, separated by distance X
% Voxel Space Dimensions
X = 200;
Y = 2; % Y boundary specularly reflective
Z = 2; % Z boundary specularly reflective

% Physical Dimensions
L = 1; % [m]: Distance between plates (Y and Z also scaled accordingly)

% Time stepping
max_itr = 200;

% Nondimensional conduction-to-radiation parameter (Modest p. 725)
N_params = [10,1,0.1,0.01,0.001]; 
N_cases = length(N_params);

% Free parameters
T1 = 600; % [K]: Temperature of plate 1 (surface at x = 1) (shouldn't matter)

% Initial nondimension temperature in medium
theta_0 = 0.75;

% Number of rays per iteration
N_rays = [1*10^5,5*10^5];%2*10^6]; 

file_name = 'ParallelPlatesCondPM';

%% Load Exact Solution
results_table = load("ParallelPlatesCondPM_Exact.mat").results_table;
x_exact = results_table{:,1};
N_param_strings = ["10","1","0.1","0.01","0.001"];
N_params_exact = double(N_param_strings);

nondim_T_exact = cell(N_cases,1);
for i = 1:N_cases
    ind_N_param = find(N_params_exact == N_params(i));
    nondim_T_exact{i} = results_table{:,strcat("N = ",N_param_strings(ind_N_param))};
end

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K^4];
tau_L = 1; % Results from modest are just for tau_L = 1.0
T2 = T1/2; % [K]: Results from modest are just for T2/T1 = 0.5
eps1 = 1; % emissivity of plate 1; results from Modest are just for black plates.
eps2 = 1; % emissivity of plate 2;

%% Derived Parameters
vx_scale = L/X; % [m/vx]: Scale of voxels
size_VS = [X+2,Y,Z]; %  +2 because each plate is 1 vx thick
PM_kappa = tau_L/X; % [1/vx]: 
PM_kappa_real = PM_kappa/vx_scale; % [1/m]:
k = N_params*4*sigma*T1^3/PM_kappa_real; % [W/(m-K)]: By definition of N (Modest p 725) (note this is a vector of k values)

%% Generate Voxel Space

VS_plate1 = false(size_VS);
VS_plate1(1,:,:) = 1;
VS_plate2 = false(size_VS);
VS_plate2(end,:,:) = 1;

VS_opaq = logical(VS_plate1 + VS_plate2); % Join the two plates into 1 opaque voxel space;
VS_opaq_eps = double(VS_opaq); % Both plates are black bodies

[VS_surf_norms, VS_surf_areas, ~] = getNormalsAndSurfaceAreas(VS_opaq,1); % Get surface normals and areas

reflective_BCs = false(2,3); % Initialize reflective BCs
reflective_BCs(:,2:3) = 1; % Y and Z boundaries are reflective 

VS_PM_kappa = zeros(size_VS); % Initialize PM voxel space
VS_PM_kappa(2:(end-1),:,:) = PM_kappa;
VS_nn = ones(size_VS); % refractive index

VS_T = zeros(size_VS); % Initialize temperature voxel space
VS_T(VS_plate1) = T1;
VS_T(VS_plate2) = T2;
VS_T(logical(VS_PM_kappa)) = theta_0*T1; % Initial temperature in medium.

% Store as voxel space object
voxel_space = VoxelSpace();
voxel_space.opaque_voxels = VS_opaq;
voxel_space.opaque_emissivities = VS_opaq_eps;
voxel_space.surface_normals = VS_surf_norms;
voxel_space.surface_areas = VS_surf_areas;
voxel_space.PM_absorption_coeffs = VS_PM_kappa;
voxel_space.refractive_indexes = VS_nn;
voxel_space.size = size_VS;
voxel_space.voxel_scale = vx_scale;
voxel_space.reflective_BCs = reflective_BCs;

%% Fixed temperatures for equilibrium solver
VS_T_fixed = false(size_VS); % Logical matrix of fixed temperatures
VS_T_fixed(VS_plate1) = 1;
VS_T_fixed(VS_plate2) = 1;

n_step = 0;
VS_T_old = VS_T;



VS_T_equil = cell(N_cases,1);
for i = 1:N_cases
    fprintf("Solving %d of %d    time = %0.3f \n",i,N_cases,toc(outer_tic))
    voxel_space.thermal_conductivity = k(i);
    [VS_T_equil{i}, ~, ~,count_itr(i)] = equilibriumCondRad(N_rays,VS_T_old,VS_T_fixed,voxel_space,"MaxIterations",max_itr);
end
total_simuation_time = toc(outer_tic)

%% Plot Solutions:
color_scheme = load(fullfile(root_folder,'src','data','ColorSchemes','PREC.mat')).color_scheme;
line_width = 1;


f = figure;
hold on;
xlim([0,1]);
ylim([0.5,1]);
set(gca,'LooseInset',get(gca,'TightInset'));
f.Position = [f.Position(1),f.Position(2),484,396];
xlabel({"Nondimensional location","$x/L$"},'Interpreter','latex')
ylabel({"Nondimensional temperature","$\theta$"},'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');
yticks([0.5,0.6,0.7,0.8,0.9,1.0])
x_vals = linspace(0,1,X+2);

labels = N_param_strings;
labels(1) = strcat("N = ",labels(1));
labels = strcat("$",labels,"$"); 
labels_x = [0.38,0.6,0.705,0.856,0.90];
labels_y = [0.726,0.726,0.726,0.728,0.77];

for i = 1:N_cases
    if i == 1
        p1_legend = plot(nan,nan,'Color',color_scheme(1),'Linewidth',line_width);
        p2_legend = plot(nan,nan,'--','Color','k','Linewidth',line_width);
    end

    nondim_T = mean(VS_T_equil{i},[2 3])/T1;
    plot(x_vals,nondim_T,'Color',color_scheme(1),'Linewidth',line_width);
    plot(x_exact,nondim_T_exact{i},'--','Color','k','Linewidth',line_width);
    text(labels_x(i),labels_y(i),labels(i),'Interpreter','latex')
end


legend({"Ray tracer with conduction","Reference solution"}, ...
    'Box','off', ...
    'Interpreter','latex');

fontsize(f,'increase');
fontsize(f,'increase');
