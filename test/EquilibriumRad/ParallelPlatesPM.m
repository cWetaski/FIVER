clear
close all
clc
%% Parallel Black Plates with Participating Medium at Equilibrium
% parallel, infinite plates (note that this is really a 1D problem)
% This reproduces the solution given in Modest section 14.3.1.2, Figure 14-3
% This is a simplified version of the full code, which does not solve to as high resolution and does not consider the whole set of Tau_L
outer_tic = tic;
%% Params
% Plates in YZ plane, separated by distance X
X = 100;
Y = 2; % Y and Z boundaries will be specularly reflecting
Z = 2;

N_rays = [5*10^4,1*10^6,5*10^6]; % Set number of rays (can be an array of rays)
N_par_workers = 6;

max_itr = 200; % Maximum number of iterations we allow
tau_L = [0.01,0.1,0.5,1,2,10];

% Temperature values are arbitrary.
T1 = 600; % [K]: Temperature of plate 1 (surface at x = 1) 
T2 = 300; % [K]: Temperature of plate 2 (surface at x = 1+X)

eps1 = 1; % emissivity of plate 1;
eps2 = 1; % emissivity of plate 2;
T_i = ((T1^4+T2^4)/2)^(1/4); % initial guess for PM temperature

scale_vx = 1; % [m/vx]: Scale of voxels

visualize = true; % Whether to visualize voxel space at the end

file_name = 'ParallelPlatesPM';

%% Start parallel pool
if ~isempty(gcp('nocreate')) % Delete existing parallel pool, if it's already running
    delete(gcp('nocreate'))
end
parpool('Threads',N_par_workers); % Threads profile is generally preferred if the code is written in such as way as to allow Threads profile.

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K];

%% Derived Parameters
size_VS = [X+2,Y,Z]; %  +2 because each plate is 1 vx thick, so opposing surface distance is X
PM_kappa = tau_L/X; % [1/vx]: Linear absorption coefficient
tau_L_strings = string(tau_L); % Used for plotting and retrieving exact solution
N_tau_L = length(tau_L);

%% Get folders
% Get current folder
cur_folder = pwd;
cur_folder = fileparts(cur_folder); % Just want the folder

% Get plots folder and project root folder
folders = regexp(cur_folder,'\','split');
for i = length(folders):-1:1 % move backward thru folders until you find VoxelRayTracer folder
    if folders(i) == "FIVER"
        root_folder = strjoin(folders,'\');
        full_file_path = fullfile(root_folder,'plots',file_name);
    else
        folders(i) = [];
    end
end

%% Load exact solution
exact_table = load(strcat(cur_folder,"\ParallelPlatesPM_Exact\ParallelPlatesPM_Exact.mat")).results_table; 
x_vals_nd_exact = exact_table{2:2:(end-1),1}; % Remove every 2nd row to make it more manageable performance-wise;
                                              % Remove 1st and last element because those values are just 1 and 0;
nondim_power_exact = cell(N_tau_L,1);
table_strings = strcat("\tau_L = ",tau_L_strings);
for i = 1:N_tau_L
    nondim_power_exact{i} = exact_table{2:2:(end-1),table_strings(i)}; % Don't grab 1st and last element since they are 1 and 0
end

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
VS_PM_kappa(2:(end-1),:,:) = PM_kappa(1);
VS_nn = ones(size_VS); % refractive index

VS_T = zeros(size_VS); % Initialize temperature voxel space
VS_T(VS_plate1) = T1;
VS_T(VS_plate2) = T2;
VS_T(logical(VS_PM_kappa)) = T_i; % Initial guess;

% Store as voxel space object
voxel_space = VoxelSpace();
voxel_space.opaque_voxels = VS_opaq;
voxel_space.opaque_emissivities = VS_opaq_eps;
voxel_space.surface_normals = VS_surf_norms;
voxel_space.surface_areas = VS_surf_areas;
voxel_space.PM_absorption_coeffs = VS_PM_kappa;
voxel_space.refractive_indexes = VS_nn;
voxel_space.size = size_VS;
voxel_space.voxel_scale = scale_vx;
voxel_space.reflective_BCs = reflective_BCs;

%% Fixed temperatures for equilibrium solver
VS_T_spec = false(size_VS);
VS_T_spec(VS_plate1) = 1;
VS_T_spec(VS_plate2) = 1;

N_temp = N_rays(end);
VS_T_equil = cell(N_tau_L,1);

if isempty(gcp('nocreate')) % Start parallel pool
    parpool;
end

%% Run Simulation
for i = 1:N_tau_L
    fprintf("Solving %d of %d    time = %0.3f \n",i,N_tau_L,toc(outer_tic))

    if tau_L(i) == 0.01
        N_rays(end) = N_temp*10; % Extra variance reduction for optically thin case
    elseif tau_L(i) == 0.1
        N_rays(end) = N_temp*4;
    else
        N_rays(end) = N_temp;
    end
    
    VS_PM_kappa(2:(end-1),:,:) = PM_kappa(i);
    voxel_space.PM_absorption_coeffs = VS_PM_kappa;
    VS_T_equil{i} = equilibriumRad(N_rays,VS_T,VS_T_spec,voxel_space,"MaxIterations",max_itr);
end
total_simulation_time = toc(outer_tic);
%% Plotting Results
dX = 1/X; % Non-dimensional distance between each voxel center
nondim_x = linspace(dX/2,1-dX/2,X);    % nondimensional x values have to be slightly adjusted since in voxel space the plate has a finite thickness
                                        % but in the exact solution X(1) corresponds to the surface of the plate.
Phi = cell(N_tau_L,1);
Phi_1D = cell(N_tau_L,1);
f = figure;
xlim([0,1]);
ylim([0,1]);
yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
xlabel({"Nondimensional location","$x/L$"},'Interpreter','latex')
ylabel({"Nondimensional emissive power";"$\Phi$"},'Interpreter','latex')
f.Position = [f.Position(1),f.Position(2),550,450];
hold on;

% load color_scheme
color_scheme = load(fullfile(root_folder,'src','data','ColorSchemes','PREC.mat')).color_scheme;
line_width = 1;

tau_L_labels = tau_L_strings;
tau_L_labels(1) = strcat("\tau_L = ",tau_L_labels(1));
x_ann = [0.02,0.08,0.08,0.08,0.08,0.08];
y_ann = [0.46,0.585,0.665,0.725,0.8,0.89];
for i = 1:N_tau_L
    Phi{i} = (VS_T_equil{i}(2:(end-1),:,:).^4-T2^4)/(T1^4-T2^4);
    Phi_1D{i} = mean(Phi{i},[2 3]); % mean at each x value
    plot(nondim_x,Phi_1D{i},'Color',color_scheme(1),'Linewidth',line_width);
    plot(x_vals_nd_exact,nondim_power_exact{i},'Color','k','LineStyle','--','Linewidth',line_width)
end
legend("Ray tracer","Reference solution", ...
    'Box','off',...
    'Interpreter','latex');

fontsize(f,'increase');

for i = 1:N_tau_L
    text(x_ann(i),y_ann(i),tau_L_labels(i));
end

fontsize(f,'increase');