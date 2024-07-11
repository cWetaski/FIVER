clear 
close all
clc

size_VS = [10,10,10];
radius = 5;
N_rays = 1e6;
power_flux = 1000; % W

voxel_space = VoxelSpace();
voxel_space.size = size_VS;
VS_opaq = false(size_VS-2);
VS_opaq = padarray(VS_opaq,[1,1,1],true,'both');
Vxyz = [1,1,1];
VS_eps = double(VS_opaq);
[VS_surf_norms, VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,1,Vxyz);
VS_nn = ones(size_VS);
voxel_space.opaque_voxels = VS_opaq;
voxel_space.opaque_emissivities = VS_eps;
voxel_space.surface_normals = VS_surf_norms;
voxel_space.surface_areas = VS_surf_areas;
voxel_space.voxel_scale = 1;
voxel_space.PM_absorption_coeffs = zeros(size_VS);
voxel_space.refractive_indexes = VS_nn;
voxel_space.reflective_BCs = zeros(2,3);
voxel_space.Vxyz = Vxyz;

VS_T = zeros(size_VS);

line_p1 = [size_VS(1)/2,size_VS(2)/2,1.01];
line_p2 = [size_VS(1)/2,size_VS(2)/2,size_VS(3)-1.01];
line_vec = (line_p2-line_p1);
line_vec = line_vec/norm(line_vec);
line_orthog = null(line_vec)';
line_orthog_vec = line_orthog(1,:);

ray_pos_generator = @(N_rays) line_p1 + line_vec.*rand(N_rays,1);
lambertian_ray_dir = @(phis,thetas) [sin(thetas).*cos(phis), sin(thetas).*sin(phis), cos(thetas)];
lambertian_ray_generator = @(N_rays) lambertian_ray_dir(2*pi*rand(N_rays,1),asin(sqrt(rand(N_rays,1))));
rotate_about_u = @(rays,thetas,u) [rays(:,1).*(cos(thetas)+u(1)^2*(1-cos(thetas)))+rays(:,2).*(u(1)*u(2)*(1-cos(thetas))-u(3)*sin(thetas))+rays(:,3).*(u(1)*u(2)*(1-cos(thetas))+u(2)*sin(thetas)), ...
                                   rays(:,1).*(u(1)*u(2)*(1-cos(thetas))+u(3)*sin(thetas))+rays(:,2).*(cos(thetas)+u(2)^2*(1-cos(thetas)))+rays(:,3).*(u(2)*u(3)*(1-cos(thetas))-u(1)*sin(thetas)), ...
                                   rays(:,1).*(u(3)*u(1)*(1-cos(thetas))-u(2)*sin(thetas))+rays(:,2).*(u(3)*u(2)*(1-cos(thetas))+u(1)*sin(thetas))+rays(:,3).*(cos(thetas)+u(3)^2*(1-cos(thetas)))];
transform_to_line = @(rays) transformCoord(rays,line_orthog_vec);
ray_dir_generator = @(N_rays) rotate_about_u(transform_to_line(lambertian_ray_generator(N_rays)),2*pi*rand(N_rays,1),line_vec);
ray_fun = @(N_rays) [ray_pos_generator(N_rays),ray_dir_generator(N_rays)];

voxel_space.addFlux(InternalFlux(power_flux,ray_fun));

phis = 2*pi*rand(N_rays,1);
thetas = asin(sqrt(rand(N_rays,1)));
rays_dir_lambertian = lambertian_ray_dir(phis,thetas);
rotation_values = 2*pi*rand(N_rays,1);
rays_dir = ray_dir_generator(N_rays);
phis2 = atan2(rays_dir(:,2),rays_dir(:,1));
histogram(phis2,100,'BinLimits',[-pi,pi],'Normalization','cdf')
figure
histogram(rays_dir_lambertian(:,1),100,'BinLimits',[-1,1],'Normalization','pdf')
hold on
histogram(rays_dir(:,3),100,'BinLimits',[-1,1],'Normalization','pdf','FaceAlpha',1)
x = linspace(-1,1,1000);
plot(x,2/pi*(1-x.^2).^0.5)
legend()
fprintf("%0.4f,   %0.4f,    %0.4f \n",sum(rays_dir(:,1)>0)/N_rays,sum(rays_dir(:,2)>0)/N_rays,sum(rays_dir(:,3)>0)/N_rays)
figure
histogram(rays_dir(:,1),100,'BinLimits',[-1,1],'Normalization','pdf','FaceAlpha',1)
hold on
figure
histogram(rays_dir(:,2),100,'BinLimits',[-1,1],'Normalization','pdf','FaceAlpha',1)
legend()
flux_power = 1000; % W

%% Start parallel pool
if ~isempty(gcp('nocreate')) % Delete existing parallel pool, if it's already running
    delete(gcp('nocreate'))
end
N_par_workers = 4;
parpool('Threads',N_par_workers); % Threads profile is generally preferred

[VS_Delta_Q,~,~] = radiativeHeatFlowsMC(N_rays,VS_T,voxel_space); % Not really sure how to test the accuracy of the line source distribution except by doing some involved view factor calculations which I really don't feel like doing!!

total_flux_abs = sum(VS_Delta_Q(:));

delete(gcp('nocreate'))




