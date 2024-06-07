clear 
close all
clc

size_VS = [10,10,10];
radius = 5;
N_rays = 1e6;

voxel_space = VoxelSpace();
voxel_space.size = size_VS;
VS_opaq = false(size_VS-2);
VS_opaq = padarray(VS_opaq,[1,1,1],true,'both');
VS_opaq(:,:,end) = false;
VS_eps = double(VS_opaq);
[VS_surf_norms, VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,1);
for ii = 2:size_VS(1)-1 % This is necessary since otherwise the surf norms near the corners can have negative dot product with ray emitted from the center, meaning the traversal function doesn't allow a collision.
    VS_surf_norms{ii,1,end-1} = [0,1,0];
    VS_surf_norms{ii,end,end-1} = [0,-1,0];
end
for jj = 2:size_VS(2)-1
    VS_surf_norms{1,jj,end-1} = [1,0,0];
    VS_surf_norms{end,jj,end-1} = [-1,0,0];
end
VS_nn = ones(size_VS);
voxel_space.opaque_voxels = VS_opaq;
voxel_space.opaque_emissivities = VS_eps;
voxel_space.surface_normals = VS_surf_norms;
voxel_space.surface_areas = VS_surf_areas;
voxel_space.voxel_scale = 1;
voxel_space.PM_absorption_coeffs = zeros(size_VS);
voxel_space.refractive_indexes = VS_nn;
voxel_space.reflective_BCs = zeros(2,3);

VS_T = zeros(size_VS);

origin_pos = [size_VS(1)/2,size_VS(2)/2,size_VS(2)/2];
ray_pos_generator = @(N_rays) repmat(origin_pos,[N_rays,1]); % point source
phi_generator = @(N_rays) 2*pi*rand(N_rays,1); % Modest eq 21.13b
thetas_generator = @(N_rays) acos(1-2*rand(N_rays,1)); % Modest eq 21.13c
ray_dir_fun = @(phis,thetas) [sin(thetas).*cos(phis),sin(thetas).*sin(phis),cos(thetas)];
ray_dir_generator = @(N_rays) ray_dir_fun(phi_generator(N_rays),thetas_generator(N_rays));
ray_fun = @(N_rays) [ray_pos_generator(N_rays),ray_dir_generator(N_rays)];

flux_power = 1000; % W
internal_flux = InternalFlux(voxel_space,flux_power,ray_fun);

[rays] = internal_flux.GenerateRays(N_rays);

%% Start parallel pool
if ~isempty(gcp('nocreate')) % Delete existing parallel pool, if it's already running
    delete(gcp('nocreate'))
end
N_par_workers = 4;
parpool('Threads',N_par_workers); % Threads profile is generally preferred

[VS_Delta_Q,~,~] = radiativeHeatFlowsMC(N_rays,VS_T,voxel_space);

total_flux_abs = sum(VS_Delta_Q(:));
predicted_flux_abs = 5/6*flux_power;
fprintf("Predicted flux: %0.2f   Actual flux: %0.2f \n",predicted_flux_abs,total_flux_abs)
delete(gcp('nocreate'))




