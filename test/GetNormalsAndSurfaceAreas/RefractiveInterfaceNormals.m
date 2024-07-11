clear
close
clc
format shortG
%% Gray enclosue geometry from Modest "Radiative Heat Transfer" 3ed example 5.4
% But with refractive interface halfway along the x-axis, checking if getting normals with refraction works correctly

vx_scale = 0.002; % [m/vx]: Scale of voxels
ns = 2; % neighbourhood size for surf normal determination

%% Params - from Modest
L = 0.40; % [m]:
h = 0.30; % [m]:
w = 0.1; % [m]: this will be repeating boundary condition though

%% Derived Parameters
size_VS = round([L/vx_scale+2,h/vx_scale+2,w/vx_scale]); % +2 to accouunt for 1vx thickness in each direction

%% Generate Voxel Space
VS_opaq = false(size_VS);
% Details from Modest example
VS_opaq(:,1,:) = 1; % S1
VS_opaq(1,:,:) = 1; % S2
VS_opaq(:,end,:) = 1; % S3
VS_opaq(end,:,:) = 1; % S4

VS_nn = ones(size_VS);
halfway_x = round(size_VS(1)/2);
VS_nn(1:halfway_x,:,:) = 2;
Vxyz = [1,1,1];

[VS_surf_norms, VS_surf_areas] = getNormalsAndSurfaceAreas(VS_opaq,ns,Vxyz,VS_nn);
VS_surf_norms_exact = cell(size_VS);
VS_surf_norms_exact(halfway_x,:,:) = {[1 0 0]}; % Refractive interface
VS_surf_norms_exact(halfway_x+1,:,:) = {[-1 0 0]}; % Refractive interface
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
if countdif ~= 2*ns*size_VS(3)*4 % the values should deviate at the four corners, 2*ns adjacent voxels at each corner.
        disp('Surface normals do not deviate from true values as expected');
else
    disp('Surface normals correct')
end
