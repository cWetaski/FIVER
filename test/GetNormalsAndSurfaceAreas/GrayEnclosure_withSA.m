clear
close
clc
format shortG
%% Gray enclosue
% All values taken from Modest, "Radiative Heat Transfer" 3ed Example 5.4
% Surface areas and normals are estimated using GetNormalsAndSurfaceAreas.m

%% Params - for ray-tracing
vx_scale = 0.02; % [m/vx]: Scale of voxels
N_rays = 10^7; % Set number of rays
N_tests = 1; % Number of times to trace N_rays
ns = 2; % neighbourhood size for surf normal determination

%% Params - from Modest (solution will be incorrect if these are modified)
L = 0.40; % [m]:
h = 0.30; % [m]:
w = 0.1; % [m]: this will be repeating boundary condition though
T_1_3 = 1000; % [K]
T_2_4 = 600; % [K]
eps_1_3 = 0.3; % emmisivity of S1 and S3
eps_2_4 = 0.8; % emissivity of S2 and S4

%% Constants
sigma = 5.670374419*10^(-8); % [W/m^2-K];

%% Derived Parameters
size_VS = round([L/vx_scale+2,h/vx_scale+2,w/vx_scale]); % +2 to accouunt for 1vx thickness in each direction

%% Generate Voxel Space
% hardcoded voxel space, being lazy;
VS_opaq = false(size_VS);
% Details from Modest example
VS_opaq(:,1,:) = 1; % S1
VS_opaq(1,:,:) = 1; % S2
VS_opaq(:,end,:) = 1; % S3
VS_opaq(end,:,:) = 1; % S4
VS_opaq_eps = double(VS_opaq)*eps_1_3; % sets e1,e3
VS_opaq_eps(1,:,:) = eps_2_4; % set e2
VS_opaq_eps(end,:,:) = eps_2_4;  % set e4
Vxyz = [1,1,1];

[VS_surf_norms, VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,ns,Vxyz);
VS_surf_norms_exact = cell(size_VS);
VS_surf_norms_exact(2:(end-1),1,:) = {[0 1 0]}; % S1
VS_surf_norms_exact(2:(end-1),end,:) = {[0 -1 0]}; % S3
VS_surf_norms_exact(1,2:(end-1),:) = {[1 0 0]}; % S2
VS_surf_norms_exact(end,2:(end-1),:) = {[-1 0 0]}; % S4

VS_surf_areas_exact = zeros(size_VS);
VS_surf_areas_exact(2:(end-1),1,:) = 1; % S1
VS_surf_areas_exact(2:(end-1),end,:) = 1; % S3
VS_surf_areas_exact(1,2:(end-1),:) = 1; % S2
VS_surf_areas_exact(end,2:(end-1),:) = 1; % S4

countdif = 0;
for i = 1:size_VS(1) 
    for j  = 1:size_VS(2)
        for k = 1:size_VS(3)
            if any(abs(VS_surf_norms_exact{i,j,k}-VS_surf_norms{i,j,k}) > 10^(-6))
                countdif = countdif+1;
            end
        end
    end
end
if countdif ~= 2*ns*size_VS(3)*4 % four corners, 2*ns adjacent voxels at each corner.
        disp('Surface normals do not deviate from true values as expected');
end

VS_T = zeros(size_VS);
VS_T(:,1,:) = T_1_3; % [K]: T1
VS_T(1,:,:) = T_2_4; % [K]: T2
VS_T(:,end,:) = T_1_3; % [K]: T3
VS_T(end,:,:) = T_2_4; % [K]: T4

reflective_BCs = false(2,3);
reflective_BCs(:,3) = 1; % set z boundary to be reflective

% Other stuff just for input to heat flow calculation
VS_PM_kappa = zeros(size_VS);
VS_nn = ones(size_VS);

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
voxel_space.Vxyz = Vxyz;


%% Analytic Solution (from Modest)
f12 = 1/4; % from modest
f14 = 1/4; % symmetry
f13 = 1/2; % summation rule;
Q1_analytic = -4230; % [W/m]: % From modest
Q2_analytic = 4230; % [W/m]:

%% Initialize Variables
Q_raytrace = zeros(N_tests,4);
VS_Q = cell(N_tests,1);
for n = 1:N_tests
    %% Generate Rays
    %fprintf("Test %d of %d \n",n,N_tests);
    VS_Q{n} = radiativeHeatFlowsMC(N_rays,VS_T,voxel_space);
end 
for n = 1:N_tests
    Q_raytrace(n,1) = sum(VS_Q{n}(2:(end-1),1,:),'all')/w; % [W/m]: Q1
    Q_raytrace(n,2) = sum(VS_Q{n}(1,:,:),'all')/w; % [W/m]: Q2
    Q_raytrace(n,3) = sum(VS_Q{n}(:,end,:),'all')/w; % [W/m]: Q3
    Q_raytrace(n,4) = sum(VS_Q{n}(end,:,:),'all')/w; % [W/m]: Q4
end

mean_Q_raytrace = mean(Q_raytrace,1);
stdev_Q_raytrace = std(Q_raytrace,1);

%% Calculate 3-sigma confidence interval
CI_lb = mean_Q_raytrace - 3*stdev_Q_raytrace/sqrt(N_tests);
CI_ub = mean_Q_raytrace + 3*stdev_Q_raytrace/sqrt(N_tests);
CI = [CI_lb; CI_ub];

Q1_error = (mean_Q_raytrace(1) - Q1_analytic)/Q1_analytic*100; % [%]: Error from expected value
Q2_error = (mean_Q_raytrace(2) - Q2_analytic)/Q2_analytic*100; % [%]:

fprintf('Q1 from raytracing: %0.0f W \n', mean_Q_raytrace(1))
fprintf('Q1 analytic: %0.0f W\n', Q1_analytic);
fprintf('Q1 error: %0.4f %% \n',Q1_error);
