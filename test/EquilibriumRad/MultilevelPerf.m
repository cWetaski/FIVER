clear
close all
clc
%% Parallel Black Plates with Participating Medium at Equilibrium
% parallel, infinite plates (note that this is really a 1D problem)
% This reproduces the solution given in Modest section 14.3.1.2, Figure 14-3
if isempty(gcp('nocreate'))
    parpool('Threads',6);
end
outer_tic = tic;
%% Params
% Plates in YZ plane, separated by distance X
X = 100;
Y = 2; % Y and Z boundaries will be specularly reflecting
Z = 2;

N_test = 20;
N_rays_multi = [5*10^4,5*10^5]; % Set number of rays (can be an array of rays)
N_rays_single1 = 2*10^5;
N_rays_single2 = 5*10^5;
max_itr = 200;
tau_L = [0.1];
tau_L_strings = string(tau_L); % Used for plotting and retrieving exact solution
N_tau_L = length(tau_L);


% Temperature values are arbitrary.
T1 = 600; % [K]: Temperature of plate 1 (surface at x = 1) 
T2 = 300; % [K]: Temperature of plate 2 (surface at x = 1+X)

eps1 = 1; % emissivity of plate 1;
eps2 = 1; % emissivity of plate 2;
T_i = ((T1+T2)/2); % initial guess for PM temperature

scale_vx = 1; % [m/vx]: Scale of voxels

visualize = true; % Whether to visualize voxel space at the end

file_name = 'ParallelPlatesPM';

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K];

%% Derived Parameters
size_VS = [X+2,Y,Z]; %  +2 because each plate is 1 vx thick, so opposing surface distance is X
PM_kappa = tau_L/X; % [1/vx]: Linear absorption coefficient

%% Get folders
% Get current folder
cur_folder = matlab.desktop.editor.getActiveFilename;
cur_folder = fileparts(cur_folder); % Just want the folder

% Get plots folder and project root folder
folders = regexp(cur_folder,'\','split');
for i = length(folders):-1:1 % move backward thru folders until you find VoxelRayTracer folder
    if folders(i) == "VoxelRayTracer"
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

[VS_surf_norms, VS_surf_areas, ~] = GetNormalsAndSurfaceAreas(VS_opaq,1); % Get surface normals and areas

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

N_temp = N_rays_multi(end);
VS_T_equil = cell(N_tau_L,1);

if isempty(gcp('nocreate')) % Start parallel pool
    parpool;
end

%% Run Simulation
dX = 1/X; % Non-dimensional distance between each voxel center
nondim_x = linspace(dX/2,1-dX/2,X);    % nondimensional x values have to be slightly adjusted since in voxel space since the plate has a finite thickness
                                        % but in the exact solution X(1) corresponds to the surface of the plate.


fprintf("Solving %d of %d    time = %0.3f \n",i,N_tau_L,toc(outer_tic))

if tau_L(i) == 0.01
    N_rays_multi(end) = N_temp*10; % Extra variance reduction for optically thin case
elseif tau_L(i) == 0.1
    N_rays_multi(end) = N_temp;
else
    N_rays_multi(end) = N_temp;
end

VS_PM_kappa(2:(end-1),:,:) = PM_kappa(i);
voxel_space.PM_absorption_coeffs = VS_PM_kappa;

for i = 1:N_test
    [VS_T_equil_multi,~,total_rays_multi(i)] = EquilibriumRad(N_rays_multi,VS_T,VS_T_spec,voxel_space,max_itr);
    [VS_T_equil_single1,~,total_rays_single1(i)] = EquilibriumRad(N_rays_single1,VS_T,VS_T_spec,voxel_space,max_itr);
    [VS_T_equil_single2,~,total_rays_single2(i)] = EquilibriumRad(N_rays_single2,VS_T,VS_T_spec,voxel_space,max_itr);

    Phi_multi =  (VS_T_equil_multi(2:(end-1),:,:).^4-T2^4)/(T1^4-T2^4);
    Phi_single1 = (VS_T_equil_single1(2:(end-1),:,:).^4-T2^4)/(T1^4-T2^4);
    Phi_single2 = (VS_T_equil_single2(2:(end-1),:,:).^4-T2^4)/(T1^4-T2^4);

    Phi_multi = mean(Phi_multi,[2 3]);
    Phi_single1 = mean(Phi_single1,[2 3]);
    Phi_single2 = mean(Phi_single2,[2 3]);

    Phi_exact = zeros(X,1);
    for j = 1:X
        x = nondim_x(j);
        j_exact = find(x_vals_nd_exact==x);
        if isempty(j_exact) % interpolate
            diffs = x_vals_nd_exact - x;
            diffs_pos = diffs;
            diffs_pos(diffs_pos<=0) = inf;
            [~,j_above] = min(diffs_pos);
            j_below = j_above - 1;
            nd_power_below = nondim_power_exact{1}(j_below);
            nd_power_above = nondim_power_exact{1}(j_above);
            Phi_exact(j) = nd_power_below + (nd_power_above-nd_power_below)/(x_vals_nd_exact(j_above)-x_vals_nd_exact(j_below))*(x-x_vals_nd_exact(j_below));
        else
            Phi_exact(j) = x_vals_nd_exact(j_exact);
        end
    end
    max_dev_multi(i) = max(abs(Phi_exact-Phi_multi)./Phi_exact)*100;
    max_dev_single1(i) = max(abs(Phi_exact-Phi_single1)./Phi_exact)*100;
    max_dev_single2(i) = max(abs(Phi_exact-Phi_single2)./Phi_exact)*100;
end
avg_max_dev_multi = mean(max_dev_multi)
avg_max_dev_single1 = mean(max_dev_single1)
avg_max_dev_single2 = mean(max_dev_single2)
avg_rays_multi = mean(total_rays_multi)
avg_rays_single1 = mean(total_rays_single1)
avg_rays_single2 = mean(total_rays_single2)
outer_toc = toc(outer_tic)