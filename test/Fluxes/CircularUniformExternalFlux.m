clear 
close all
clc

size_VS = [10,10,10];
radius = 5;
N_rays = 5e3;
center_loc = [size_VS(1)/2, size_VS(2)/2,0];
circle_loc = @(theta,r) [cos(theta),sin(theta),0*theta].*r + center_loc;
ray_pos_inv_CDF = @(N_rays) circle_loc(rand(N_rays,1)*2*pi,sqrt(rand(N_rays,1))*radius);
ray_dir_inv_CDF = @(N_rays) repmat([0,0,1],[N_rays,1]);
ray_fun = @(N_rays) [ray_pos_inv_CDF(N_rays),ray_dir_inv_CDF(N_rays)];
voxel_space = VoxelSpace();
voxel_space.size = size_VS;

VS_opaq = false(size_VS);
VS_opaq(:,end,:) = true;
VS_eps = double(VS_opaq);
[VS_surf_norms, VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,1);
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


voxel_space.addFlux(ExternalFlux(1000,ray_fun,'bottom',size_VS));


[rays_pos, rays_dir] = voxel_space.fluxes{end}.GenerateRays(N_rays);
scatter(rays_pos(:,1),rays_pos(:,3),'Marker','.')

%% Start parallel pool
if ~isempty(gcp('nocreate')) % Delete existing parallel pool, if it's already running
    delete(gcp('nocreate'))
end
N_par_workers = 4;
parpool('Threads',N_par_workers); % Threads profile is generally preferred if the code is written in such as way as to allow Threads profile.

[VS_Delta_Q,~,~] = radiativeHeatFlowsMC(N_rays,VS_T,voxel_space);
sum(VS_Delta_Q(:))

delete(gcp('nocreate'))


