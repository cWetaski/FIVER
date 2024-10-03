% ©2024 ETH Zurich, Charles Wetaski, Sebastian Sas Brunser, Emiliano Casati
%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-09-06

function [VS_surf_norms,VS_surf_areas,inds] = getNormalsAndSurfaceAreas(VS_opaq,vx_scale,ns,VS_nn)
    %GETNORMALSANDSURFACEAREAS Returns the normal vector and surface area (in vx^2) of each external surface voxel in VS
    %   INPUTS:
    %       VS_opaq (3D logical):           Voxel space of opaque voxels
    %       ns (scalar double (int)):       Neighbourhood size for normal estimation   
    %       VS_nn (3D double (nn >= 1)):    (OPTIONAL) Voxel space of refractive indices.
    %                                           For problems with refraction interfaces.
    %       Vxyz (1x3 double)               Relative size of voxel in dimension. 
    %   OUTPUTS:
    %       VS_surf_norms (3D double):      Matrix of (normalized) surface normals of each surface voxels
    %       VS_surf_areas (3D cell) [vx^2]: Vector of surface areas of each surface voxel 
    %       inds (Nx1 double):              Vector of linear indices of each surface voxel in the voxel space
    %
    size_VS = size(VS_opaq);
    VS_surf_norms = cell(size_VS);
    VS_surf_areas = zeros(size_VS);
    A_norm = [vx_scale(2)*vx_scale(3),vx_scale(1)*vx_scale(3),vx_scale(1)*vx_scale(2)]; % Area of planes normal to each direction [vx^2]

    if nargin == 4 % Want normals at refractive interfaces as well
        nn_vals = unique(VS_nn); % Get unique refractive index values
        N_nn = length(nn_vals); % Get number of unique refractive index values
        if N_nn > 1 % Need more than 1 value of nn to have refractive interfaces
            for i = 1:N_nn % iterate over different refractive indexes
                VS_region_i = (VS_nn == nn_vals(i)); % get current region
                VS_surf_norms_i = getNormalsAndSurfaceAreas(VS_region_i,vx_scale,ns); % Recursively call surface norms on refraction regions
                VS_surf_norms(VS_region_i) = VS_surf_norms_i(VS_region_i); % Assign to surf norms cell array 
            end
        end
    end
    % Now get surface normals for opaque voxels
    VS_surf_external = countConnectivity(VS_opaq); % (3D double (int)): get the 6-connectivity of each voxel (boundaries are considered solid)
    zero_connected_voxels = (VS_surf_external==0) & VS_opaq;
    N_zero_connected = sum(zero_connected_voxels(:));
    if N_zero_connected > 0
        warning("Domain contains %d zero-connected voxels, this may lead to issues",N_zero_connected);
        VS_surf_external(zero_connected_voxels)=1; % We still try to obtain a surface normal for zero-connected voxels
    end
    VS_surf_external(VS_surf_external == 6) = 0; % remove 6-connected (i.e., internal) voxels -> external surface voxel remain

    [inds] = find(VS_surf_external); % (1D col double (int)): get linear indices corresponding to external surface voxels
    N_surf_voxels = length(inds); % (scalar double (int)): Count number of external surface voxels
    
    for i = 1:N_surf_voxels
        cur_norm = getNormalVector(inds(i),VS_opaq,ns,vx_scale); % estimate the normal vector (vector is normalized)
        if any(isnan(cur_norm))
            [x,y,z] = ind2sub(size_VS,inds(i));     
            error("nan normal vector at [%d %d %d]: you can try increasing the neighbourhood size or you can double the resolution with repelem(VS_opaq,[2,2,2])",x,y,z)
        end
        
        for j = 1:3
            if abs(cur_norm(j)) > 0 && abs(cur_norm(j)) < 10^(-9) % Just a rounding error
                cur_norm(j) = 0;
                cur_norm = cur_norm/norm(cur_norm);
            end
        end
        VS_surf_norms{inds(i)} = cur_norm; % store normal vector
        VS_surf_areas(inds(i)) = 1/max(abs(cur_norm)./A_norm); % estimate the area (source: Flin et al. 2005, "Adaptive Estimation of Normals and Surface Area ...")
    end
end