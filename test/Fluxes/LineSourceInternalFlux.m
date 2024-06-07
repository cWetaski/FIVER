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
%VS_opaq(:,:,end) = false;
VS_eps = double(VS_opaq);
[VS_surf_norms, VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,1);
% for ii = 2:size_VS(1)-1 % This is necessary since otherwise the surf norms near the corners can have negative dot product with ray emitted from the center, meaning the traversal function doesn't allow a collision.
%     VS_surf_norms{ii,1,end-1} = [0,1,0];
%     VS_surf_norms{ii,end,end-1} = [0,-1,0];
% end
% for jj = 2:size_VS(2)-1
%     VS_surf_norms{1,jj,end-1} = [1,0,0];
%     VS_surf_norms{end,jj,end-1} = [-1,0,0];
% end
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
ray_dir_generator = @(N_rays) rotate_about_u(transform_to_line(lambertian_ray_generator(N_rays)),2*pi*rand(N_rays,1),line_orthog_vec);
ray_fun = @(N_rays) [ray_pos_generator(N_rays),ray_dir_generator(N_rays)];

phis = 2*pi*rand(N_rays,1);
thetas = asin(sqrt(rand(N_rays,1)));
rays_dir_lambertian = lambertian_ray_dir(phis,thetas);
rotation_values = 2*pi*rand(N_rays,1);
rays_dir_transformed = transform_to_line(rays_dir_lambertian);
rays_dir = rotate_about_u(rays_dir_transformed,rotation_values,line_vec);
phis2 = atan2(rays_dir(:,2),rays_dir(:,1));
histogram(phis2,100,'BinLimits',[-pi,pi],'Normalization','cdf')
figure
histogram(rays_dir_lambertian(:,1),100,'BinLimits',[-1,1],'Normalization','pdf')
hold on
histogram(rays_dir(:,3),100,'BinLimits',[-1,1],'Normalization','pdf','FaceAlpha',1)
x = linspace(-1,1,1000);
plot(x,2/pi*(1-x.^2).^0.5)
legend()

%hold on
%histogram(thetas2,100,'BinLimits',[0,pi/2],'Normalization','cdf')
fprintf("%0.4f,   %0.4f,    %0.4f \n",sum(rays_dir(:,1)>0)/N_rays,sum(rays_dir(:,2)>0)/N_rays,sum(rays_dir(:,3)>0)/N_rays)

figure
histogram(cos(phis),100,'BinLimits',[-1,1],'Normalization','pdf','FaceAlpha',1)
hold on
f = @(z) 1/pi*sqrt(1-z.^2);
v = linspace(-1,1,1000);
plot(v,1./(pi*sqrt(1-v.^2)));
legend()
figure
histogram(sin(thetas),100,'BinLimits',[0,1],'Normalization','pdf','FaceAlpha',1);
hold on
u = linspace(0,1,1000);
plot(u,2*u);
flux_power = 1000; % W
internal_flux = InternalFlux(voxel_space,flux_power,ray_fun); % Returns the internal flux but also assigns it to the voxel space

[rays] = internal_flux.GenerateRays(N_rays);

%% Start parallel pool
if ~isempty(gcp('nocreate')) % Delete existing parallel pool, if it's already running
    delete(gcp('nocreate'))
end
N_par_workers = 4;
parpool('Threads',N_par_workers); % Threads profile is generally preferred

[VS_Delta_Q,~,~] = radiativeHeatFlowsMC(N_rays,VS_T,voxel_space);

total_flux_abs = sum(VS_Delta_Q(:))

delete(gcp('nocreate'))




