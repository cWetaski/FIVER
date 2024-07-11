%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07

function [VS_surf_norms,VS_surf_areas,inds] = getNormalsAndSurfaceAreas(VS_opaq,ns,Vxyz,VS_nn)
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
    inds_sparse = zeros(size_VS);
    A_norm = [Vxyz(2)*Vxyz(3),Vxyz(1)*Vxyz(3),Vxyz(1)*Vxyz(2)]; % Area of planes normal to each direction [vx^2]

    if nargin == 4 % Want normals at refractive interfaces as well
        nn_vals = unique(VS_nn); % Get unique refractive index values
        N_nn = length(nn_vals); % Get number of unique refractive index values

        for i = 1:N_nn % iterate over different refractive indexes
            VS_region_i = (VS_nn == nn_vals(i)); % get current region
            VS_surf_norms_i = getNormalsAndSurfaceAreas(VS_region_i,ns); % Recursively call surface norms on refraction regions
            VS_surf_norms(VS_region_i) = VS_surf_norms_i(VS_region_i); % Assign to surf norms cell array 
        end
    end
    % Now get surface normals for opaque voxels
    VS_surf_external = countConnectivity(VS_opaq); % (3D double (int)): get the 6-connectivity of each voxel (boundaries are considered solid)
    VS_surf_external(VS_surf_external == 6) = 0; % remove 6-connected (i.e., internal) voxels -> external surface voxel remain

    inds = find(VS_surf_external); % (1D col double (int)): get linear indices corresponding to external surface voxels
    
    N_surf_voxels = length(inds); % (scalar double (int)): Count number of external surface voxels
    
    surf_norms_lin = zeros(N_surf_voxels,3);
    for i = 1:N_surf_voxels
        cur_norm = getNormalVector(inds(i),VS_opaq,ns,Vxyz); % estimate the normal vector (vector is normalized)
        for j = 1:3
            if abs(cur_norm(j)) > 0 && abs(cur_norm(j)) < 10^(-9) % Just a rounding error
                cur_norm(j) = 0;
                cur_norm = cur_norm/norm(cur_norm);
                break;
            end
        end
        VS_surf_norms{inds(i)} = cur_norm; % store normal vector
        inds_sparse(inds(i)) = i;
        surf_norms_lin(i,:) = cur_norm;
        VS_surf_areas(inds(i)) = 1/max(abs(cur_norm)./A_norm); % estimate the area (source: Flin et al. 2005, "Adaptive Estimation of Normals and Surface Area ...")
    end
    inds_sparse = sparse(inds_sparse(:));
end