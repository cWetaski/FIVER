%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07

function [ray_pos,ray_dir] = generateSurfaceRays(pos_vx,surf_norm,N_rays)
    % GENERATESURFACERAY Generates a ray from a surface voxel based on the voxel coordinate and surface normal 
    % INPUT:
    %   pos_vx (1x3 double (int)):      Voxel coordinate in voxel space
    %   surf_norm (1x3 double):         Local normal vector of the surface
    %   N_rays (scalar double (int)):   (OPTIONAL) Number of rays to generate (default = 1)
    % OUTPUT:
    %   ray_pos (Nx3 double):           Ray position (N = N_rays)
    %   ray_dir (Nx3 double):           Ray direction (N = N_rays)
    %
    if nargin == 2
        N_rays = 1;
    end
    projection_distance_factor = 1-1e-10; % This factor determines how much past the border of the voxel the ray will originate
                                       % from, setting it equal to 1 means that the ray will originate on the border of the
                                       % voxel, however, this can lead to incorrect solutions do to self intersection when
                                       % dealing with curved surface. However, setting it greater than 1 will lead to
                                       % incorrect results for PM voxels near surfaces since a portion of the ray traversal
                                       % will always be skipped. Keeping it set to 0.99 is safest (the ray is emitted from just within the voxel boundary)
    % Generate ray position
    proj_dist = 1/2*1/max(abs(surf_norm)); % Distance to border of voxel from center of voxel.
    %orthovecs = null(surf_norm(:).')';
    ray_pos = (pos_vx - 0.5)+ proj_dist*surf_norm*projection_distance_factor;
    ray_pos = repmat(ray_pos,[N_rays,1]); % repeat N copies
    %ray_pos = ray_pos + 0.5*(0.5-rand(N_rays,1))*orthovecs(1,:) + 0.5*(0.5-rand(N_rays,1))*orthovecs(2,:);
    % Generate emission direction (lambertian emitter)
    thetas = asin(sqrt(rand(N_rays,1)));    % Modest Eq. 8.42
    phis=2*pi*rand(N_rays,1);               % Modest Eq. 8.41
    dir_local = [sin(thetas).*cos(phis), sin(thetas).*sin(phis), cos(thetas)]; % Direction of random vector with the normal/local coordinates
    ray_dir = transformCoord(dir_local,surf_norm); % Transform directions to global
end

