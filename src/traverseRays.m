% ©2024 ETH Zurich, Charles Wetaski, Sebastian Sas Brunser, Emiliano Casati
%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-09-06

function [rays_end_pos,rays_events] = traverseRays(rays_pos_start,rays_dir_start,VS_opaq,VS_opaq_eps,VS_surf_norms,VS_PM_kappa,VS_nn,reflective_BCs,vx_scale,size_VS)
    % TRAVERSERAYS - Traces a set of rays from starting position and direction through voxel space until absorption or
    %                 exiting the vox\el space. 
    %                   
    % INPUTS:
    %   rays_pos_start (Nx3 double) [vx]:               Positions where rays originate
    %   rays_dir_start (Nx3 double) [vx]:               Directions of rays
    %   VS_opaq (3D logical):                           Stores whether a voxel is opaque (i.e., not PM and not empty)     
    %   VS_opaq_eps (3D double (in [0,1])) [-]:         Stores emissivity of opaque voxels
    %   VS_surf_norms (3D cell):                        Stores 1x3 vectors of surface normals. Empty cell if not a
    %                                                       surface voxel.
    %   VS_PM_kappa (3D double (>= 0)) [1/vx]:          Stores the linear absorption coefficient of PM voxels.
    %   VS_nn   (3D double (>= 1)) [1/vx]:              Stores refractive index at each location
    %   reflective_BCs (2x3 logical):                   Boundary conds: Rows are lower/upper bound, cols are XYZ. A
    %                                                       reflective boundary reflects the ray specularly.
    %   vx_scale (1x3 double) [m/vx]                    Size of voxels in each dimension
    %   size_VS (1x3 double (int, sz >= 1)) [vx]:       Size of voxel space in each dimension, although this can easily
    %                                                       be determined in the function, it is ~5% faster to pass it.
    % OUTPUTS:
    %   rays_end_pos (Nx3 double (int)):                End positions. For exit events, position will correspond to voxel 
    %                                                   coordinate outside of Voxel Space.
    %   event (double (int)):                           Different possible events (listed below) use numerical encoding
    %                                                   since it is much faster when manipulating arrays of events.
    %       1 == exit: Exited the voxel space
    %       2 == medium absorption: Absorbed by medium 
    %       3 == surface absorption: Absorbed by surface
    %       4 == medium self-absorption: Absorbed by medium from which it was emitted without ever exiting
    %% Pre-preamble
    acceptance_angle = 0; % A surface collision only occurs if the dot product of the surface normal vector and the ray direction vector is less than this value (default == 0).
    max_loops = 100; % Set maximum number of loops (each iteration of the outer loop indicates a new reflection or refraction event)
    N_rays = size(rays_pos_start,1);
    rays_end_pos = zeros(N_rays,3);
    rays_events = zeros(N_rays,1);
   %% Start parallel pool
    if isempty(gcp('nocreate')) % If there's no parpool already created
        parpool('Threads'); % Threads profile is generally preferred if the code is written in such as way as to allow Threads profile.
    %else
    %    delete(gcp) % No parallel!
    end
    parfor nn_ray = 1:N_rays
        %% Initialize temporary variables to avoid warnings:
        XYZ = 0; end_pos = [-1,-1,-1]; event = 0; tau_ray = 0; tau_acc = 0; dist_frac_b = 0; dist_frac_c = 0; last_inc = 0; nn2 = 0; i_min = 0;
        
        %% Preamble Part 1
        rtrn = false;
        ray_dir = rays_dir_start(nn_ray,:);
        ray_pos = rays_pos_start(nn_ray,:);
        sgn = sign(ray_dir);
        if any(ray_pos == floor(ray_pos)) % Check if starting on boundary
            % eps == machine epsilon
            ray_pos = ray_pos + eps*sgn.*ray_pos; % Give small bump to get off boundary (scale bump with position for floating point reasons (i.e., floor(1-eps) == 0, but floor(10-eps) == 10)
        end 
        % Check starting position is inside voxel space
        for ii = 1:3 
            if ray_pos(ii) > size_VS(ii) % outside maximum boundary
                if reflective_BCs(2,ii) && ray_pos(ii) - size_VS(ii) < 1
                                                                % if only a small distance over reflective boundary
                   ray_pos(ii) = 2*size_VS(ii) - ray_pos(ii);     % reflect over boundary
                   ray_dir(ii) = -ray_dir(ii);
                else
                    end_pos = floor(ray_pos);
                    event = 1; % exit
                    rtrn = true;
                end
            elseif ray_pos(ii) < 0 % outside minimum boundary
                if reflective_BCs(1,ii) && all(ray_pos > -1)
                                                            % if only a small distance over reflective boundary
                    ray_pos(ii) = -ray_pos(ii);           % reflect over boundary
                    ray_dir(ii) = -ray_dir(ii);
                else
                    end_pos = floor(ray_pos);
                    event = 1; % exit ;
                    rtrn = true;
                end
            end
        end
        
        first_pass_bool = true;
        loop_count = 0;
        XYZ_0 =  min(floor(ray_pos)+1,size_VS); % accounts for starting coordinate which is on voxel space boundary;
        
        while ~rtrn && loop_count < max_loops % While loop since the preamble is repeated if a ray is reflected diffusely or refracted
            loop_count = loop_count+1;
            %% Preamble part 2 - This part is repeated when there is a refraction or a diffuse surface reflection 
            % Determine signs of direction vector
            % if reflection_refraction_count >= 100 
            %     fprintf('Ray Stuck in Traversal after %d reflections/refractions - Returning \n',reflection_refraction_count)
            %     end_pos = XYZ;
            %     event = 2; % Just assume its in PM
            %     % VisualizeVoxelSpace(VS_nn-1);
            %     % hold on
            %     % 
            %     % 
            %     % lb = ray_pos - 3;
            %     % ub = ray_pos + 3;
            %     % xlim([lb(1) ub(1)])
            %     % ylim([lb(2) ub(2)])
            %     % zlim([lb(3) ub(3)])
            %     % scatter3(ray_pos(1),ray_pos(2),ray_pos(3),50,'g','filled')
            %     % quiver3(ray_pos(1),ray_pos(2),ray_pos(3),2*ray_dir(1),2*ray_dir(2),2*ray_dir(3))
            %     return
            % end
        
            % Get starting distance from starting voxel boundary;
            sgn = sign(ray_dir);
            xs = (sgn.*(rem(ray_pos,1)-0.5)+0.5); % get remainder from boundary of starting voxel (note this fails if sgn = 0, but that should never occur for a real simulation.
            % Using a trick trick so that if sgn = 1: xs = rem, but if s = -1: xs = 1 - rem, since the direction is switched
        
            % Determine driving axis
            abs_ray_dir = abs(ray_dir);
            % Note this does not sort the axes; the secondary and ternary axes form a right handed coordinate system with the driving axis
            % This is faster than using MATLAB's built in sort or max functions (in my tests the built-in functions took 
            % roughly twice as long.
            if abs_ray_dir(1)*vx_scale(2) > vx_scale(1)*abs_ray_dir(2)
                if abs_ray_dir(1)*vx_scale(3)  > vx_scale(1)*abs_ray_dir(3) % x is driving axis
                    ia = 1;
                    ib = 2;
                    ic = 3;
                else % z is driving axis
                    ia = 3;
                    ib = 1;
                    ic = 2;
                end
            elseif abs_ray_dir(2)*vx_scale(3) > abs_ray_dir(3)*vx_scale(2)
                ia = 2;
                ib = 3;
                ic = 1;
            else
                ia = 3;
                ib = 1;
                ic = 2;
            end
        
            % Get default distance (i.e., the distance the ray traverses per increment of the driving axis)
            ray_dir_scaled = abs_ray_dir/abs_ray_dir(ia); % Divide by magnitude of driving axis direction da
            default_dist = norm(ray_dir_scaled)*vx_scale(ia);
            
            % Get increments of secondary axis and ternary axis (per driving axis increment)
            Db = ray_dir_scaled(ib)*vx_scale(ia)/vx_scale(ib); % (db/da) in paper this is Vx*dy/dx, -
            Dc = ray_dir_scaled(ic)*vx_scale(ia)/vx_scale(ic); % (dc/da)
        
            pb = xs(ib) - 1 + (1-xs(ia))*Db; % this corresponds to p_xy' in Liu et al paper (for x driving axis)
            pc = xs(ic) - 1 + (1-xs(ia))*Dc; % (p_xz')
        
            if first_pass_bool % Don't want to reset this if ray is diffusely reflected or refracted
                XYZ = XYZ_0; % On the first pass, we compute this before entering the while-loop
                % Generate optical distance that the ray can traverse before being absorbed
                tau_ray = -log(rand); % Modest, Eq. 21.19
                % Initialize optical depth accumulator
                tau_acc = -xs(ia)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); % subtraction to account for starting partway thru the voxel
                first_pass_bool = false; % We will not repeat this for a reflected/refracted ray+
            else % Not first pass, have to determine XYZ again.
                % XYZ stores the current voxel coordinate position (i.e., XYZ(1) = X coordinate, XYZ(2) = Y coordinate, etc.)
                XYZ =  min(floor(ray_pos)+1,size_VS); % accounts for starting coordinate which is on voxel space boundary;
            end
            
            % Get current medium 
            nn1 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
            
            surface_collision_bool = false;
            interior_collision_bool = false;
        
            %% RAY TRAVERSAL    
            while true % Traversal loop (breaks upon collision with opaque surface or medium interface)
                if pb > 0   % Hierarchy 1 
                    %% increment XYZ(ib) before XYZ(ia)
                    dist_frac_b = pb/Db;        % compute fraction of distance travelled after crossing XYZ(ib) border
                    if pc > 0 % Hierarchy 2
                        %% increment XYZ(ic) before XYZ(ia)
                        dist_frac_c = pc/Dc;        % compute fraction of distance travelled after crossing XYZ(ic) border
                        if dist_frac_b > dist_frac_c % Hierarchy 3
                            %% increment XYZ(ib) before XYZ(ic)
                                                        % Increment tau up to crossing XYZ(ib) boundary
                            tau_acc = tau_acc + (1-dist_frac_b)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                            if tau_acc > tau_ray        
                                                        % absorbed before crossing boundary
                                event = 2; % medium absorption
                                end_pos = XYZ;
                                rtrn = true;
                                break
                            end
        
                            % Increment XYZ(ib)
                            XYZ(ib) = XYZ(ib) + sgn(ib);
                            pb = pb - 1; % Reset pb
        
                            % Check exit or collision
                            if XYZ(ib) > size_VS(ib) 
                                                        % exited at maximum
                                if reflective_BCs(2,ib) 
                                    sgn(ib) = -sgn(ib);         % flip direction
                                    ray_dir(ib) = -ray_dir(ib);
                                    XYZ(ib) = size_VS(ib);    % undo previous increment
                                else
                                    end_pos = XYZ;            % Note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif XYZ(ib) < 1 
                                                        % exited at minimum
                                if reflective_BCs(1,ib)     
                                    sgn(ib) = -sgn(ib);         % flip direction
                                    ray_dir(ib) = -ray_dir(ib);
                                    XYZ(ib) = 1;              % undo increment
                                else
                                    end_pos = XYZ;            % note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                                surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                                if isempty(surf_norm) 
                                    last_inc = ib;
                                    interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                                    break
                                elseif dot(surf_norm,ray_dir) < acceptance_angle
                                    surface_collision_bool = true;   % collision with opaque surface voxel
                                    break
                                end
                            else
                                % Check refraction
                                nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                                if nn2 ~= nn1 % if refractive index changes
                                    surface_collision_bool = false; % Medium refractive interface
                                    last_inc = ib; % Last increment was ib
                                    break
                                end
                            end
                                
                            % Increment tau up to crossing XYZ(ic) boundary
                            tau_acc = tau_acc + (dist_frac_b - dist_frac_c)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3));
                            if tau_acc > tau_ray
                                                            % absorbed before crossing boundary
                                event = 2; % medium absorption
                                end_pos = XYZ;
                                rtrn = true;
                                break
                            end
        
                            % Increment XYZ(ic)
                            XYZ(ic) = XYZ(ic) + sgn(ic);
                            pc = pc - 1; % Reset pc
        
                            % Check exit or collision
                            if XYZ(ic) > size_VS(ic) 
                                                        % exited at maximum
                                if reflective_BCs(2,ic)     
                                    sgn(ic) = -sgn(ic);         % flip direction
                                    ray_dir(ic) = -ray_dir(ic);
                                    XYZ(ic) = size_VS(ic);    % undo increment
                                else
                                    end_pos = XYZ;            % note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif XYZ(ic) < 1 
                                                        % exited at minimum
                                if reflective_BCs(1,ic)     
                                    sgn(ic) = -sgn(ic);         % flip direction
                                    ray_dir(ic) = -ray_dir(ic);
                                    XYZ(ic) = 1;              % undo increment
                                else
                                    end_pos = XYZ;            % note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                                surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                                if isempty(surf_norm) 
                                    last_inc = ic;
                                    interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                                    break
                                elseif dot(surf_norm,ray_dir) < acceptance_angle
                                    surface_collision_bool = true;   % collision with opaque surface voxel
                                    break
                                end
                            else
                                % Check refraction
                                nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                                if nn2 ~= nn1 % if refractive index changes
                                    surface_collision_bool = false; % Medium refractive interface
                                    last_inc = ic; % Last increment was ic
                                    break
                                end
                            end
        
                            % Increment tau up to crossing XYZ(ia) booundary
                            tau_acc = tau_acc + dist_frac_c*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3));
                            if tau_acc > tau_ray
                                                        % absorbed before crossing boundary
                                event = 2; % medium absorption
                                end_pos = XYZ;
                                rtrn = true;
                                break
                            end
        
                        else % Hierarchy 3
                            %% increment XYZ(ic) before XYZ(ib)
                                                        % Increment tau up to crossing XYZ(ic) boundary
                            tau_acc = tau_acc + (1-dist_frac_c)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                            if tau_acc > tau_ray 
                                                        % absorbed before crossing boundary
                                event = 2; % medium absorption
                                end_pos = XYZ;
                                rtrn = true;
                                break;
                            end
        
                            % Increment XYZ(ic)
                            XYZ(ic) = XYZ(ic) + sgn(ic); 
                            % Reset pc
                            pc = pc - 1;     
        
                            % Check exit or collision
                            if XYZ(ic) > size_VS(ic) 
                                                        % exited at maximum
                                if reflective_BCs(2,ic)     
                                    sgn(ic) = -sgn(ic);         % flip direction
                                    ray_dir(ic) = -ray_dir(ic);
                                    XYZ(ic) = size_VS(ic);    % undo increment
                                else
                                    end_pos = XYZ;            % note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break;
                                end
                            elseif XYZ(ic) < 1 
                                                        % exited at minimum
                                if reflective_BCs(1,ic)     
                                    sgn(ic) = -sgn(ic);         % flip direction
                                    ray_dir(ic) = -ray_dir(ic);
                                    XYZ(ic) = 1;              % undo increment
                                else
                                    end_pos = XYZ;            % note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                                surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                                if isempty(surf_norm) 
                                    last_inc = ic;
                                    interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                                    break
                                elseif dot(surf_norm,ray_dir) < acceptance_angle
                                    surface_collision_bool = true;   % collision with opaque surface voxel
                                    break
                                end
                            else
                                % Check refraction
                                nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                                if nn2 ~= nn1 % if refractive index changes
                                    surface_collision_bool = false; % Medium refractive interface
                                    last_inc = ic; % Last increment was ic
                                    break
                                end
                            end
        
                            % Increment tau up to crossing XYZ(ib) boundary
                            tau_acc = tau_acc + (dist_frac_c-dist_frac_b)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                            if tau_acc > tau_ray 
                                                        % absorbed before crossing boundary
                                event = 2; % medium absorption
                                end_pos = XYZ;
                                rtrn = true;
                                break
                            end
                            
                            % Increment XYZ(ib)
                            XYZ(ib) = XYZ(ib) + sgn(ib);  
                            pb = pb - 1; % Reset pb
                            
                            % Check exit or collision
                            if XYZ(ib) > size_VS(ib) 
                                                        % exited at maximum
                                if reflective_BCs(2,ib) 
                                    sgn(ib) = -sgn(ib);         % flip direction
                                    ray_dir(ib) = -ray_dir(ib);
                                    XYZ(ib) = size_VS(ib);    % undo previous increment
                                else
                                    end_pos = XYZ;            % Note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif XYZ(ib) < 1 
                                                        % exited at minimum
                                if reflective_BCs(1,ib)     
                                    sgn(ib) = -sgn(ib);         % flip direction
                                    ray_dir(ib) = -ray_dir(ib);
                                    XYZ(ib) = 1;              % undo increment
                                else
                                    end_pos = XYZ;            % note that X will outside of voxel_space
                                    event = 1; % exit
                                    rtrn = true;
                                    break
                                end
                            elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                                surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                                if isempty(surf_norm) 
                                    last_inc = ib;
                                    interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                                    break
                                elseif dot(surf_norm,ray_dir) < acceptance_angle
                                    surface_collision_bool = true;   % collision with opaque surface voxel
                                    break
                                end
                            else
                                % Check refraction
                                nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                                if nn2 ~= nn1 % if refractive index changes
                                    surface_collision_bool = false; % Medium refractive interface
                                    last_inc = ib; % Last increment was ib
                                    break
                                end
                            end
        
                            % Increment tau up to crossing XYZ(ia) boundary
                            tau_acc = tau_acc + dist_frac_b*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                            if tau_acc > tau_ray 
                                                        % absorbed before crossing boundary
                                event = 2; % medium absorption
                                end_pos = XYZ;
                                rtrn = true;
                                break
                            end
                        end % Hierarchy 3
                                 
                    else % Hierarchy 2
                        %% Increment only XYZ(ib) before XYZ(ia)
                        
                        % Increment tau up to crossing XYZ(ib) boundary
                        tau_acc = tau_acc + (1-dist_frac_b)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                        if tau_acc > tau_ray 
                                                    % absorbed before crossing boundary
                            event = 2; % medium absorption
                            end_pos = XYZ;
                            rtrn = true;
                            break
                        end
        
                        % Increment XYZ(ib)
                        XYZ(ib) = XYZ(ib) + sgn(ib);  
                        pb = pb - 1; % Reset pb
                        
                        % Check exit or reflection
                        if XYZ(ib) > size_VS(ib) 
                                                    % exited at maximum
                            if reflective_BCs(2,ib) 
                                sgn(ib) = -sgn(ib);         % flip direction
                                ray_dir(ib) = -ray_dir(ib);
                                XYZ(ib) = size_VS(ib);    % undo previous increment
                            else
                                end_pos = XYZ;            % Note that X will be outside of voxel_space
                                event = 1; % exit
                                rtrn = true;
                                break
                            end
                        elseif XYZ(ib) < 1 
                                                    % exited at minimum
                            if reflective_BCs(1,ib)     
                                sgn(ib) = -sgn(ib);         % flip direction
                                ray_dir(ib) = -ray_dir(ib);
                                XYZ(ib) = 1;              % undo increment
                            else
                                end_pos = XYZ;            % note that X will be outside of voxel_space
                                event = 1; % exit
                                rtrn = true;
                                break
                            end
                        elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                            surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                            if isempty(surf_norm) 
                                last_inc = ib;
                                interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                                break
                            elseif dot(surf_norm,ray_dir) < acceptance_angle
                                surface_collision_bool = true;   % collision with opaque surface voxel
                                break
                            end
                        else
                            % Check refraction
                            nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                            if nn2 ~= nn1 % if refractive index changes
                                surface_collision_bool = false; % Medium refractive interface
                                last_inc = ib; % Last increment was ib
                                break
                            end
                        end
                        
                        % Increment tau up to crossing XYZ(ib) boundary
                        tau_acc = tau_acc + (dist_frac_b)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                        if tau_acc > tau_ray 
                                                    % absorbed before crossing boundary
                            event = 2; % medium absorption
                            end_pos = XYZ;
                            rtrn = true;
                            break
                        end
                    end % Hierarchy 2               
                else % Hierarchy 1
                    %% Do not increment XYZ(ib)
                    if pc > 0
                        %% increment only XYZ(ic) before XYZ(ia)
        
                        % compute fraction of distance travelled after crossing XYZ(ic) boundary
                        dist_frac_c = pc/Dc;        
                                                    
                        % Increment tau up to crossing XYZ(ic) boundary
                        tau_acc = tau_acc + (1-dist_frac_c)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                        if tau_acc > tau_ray 
                            % absorbed before crossing boundary
                            event = 2; % medium absorption
                            end_pos = XYZ;
                            rtrn = true;
                            break
                        end
        
                        % Increment XYZ(ic)
                        XYZ(ic) = XYZ(ic) + sgn(ic);
                        pc = pc - 1; % Reset ic
        
                        % Check exit or collision
                        if XYZ(ic) > size_VS(ic) 
                                                    % exited at maximum
                            if reflective_BCs(2,ic)     
                                sgn(ic) = -sgn(ic);         % flip direction
                                ray_dir(ic) = -ray_dir(ic);
                                XYZ(ic) = size_VS(ic);    % undo increment
                            else
                                end_pos = XYZ;            % not0e that X will outside of voxel_space
                                event = 1; % exit
                                rtrn = true;
                                break
                            end
                        elseif XYZ(ic) < 1 
                                                    % exited at minimum
                            if reflective_BCs(1,ic)     
                                sgn(ic) = -sgn(ic);         % flip direction
                                ray_dir(ic) = -ray_dir(ic);
                                XYZ(ic) = 1;    % undo increment
                            else
                                end_pos = XYZ;            % note that X will outside of voxel_space
                                event = 1; % exit
                                rtrn = true;
                                break
                            end
                            elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                                surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                                if isempty(surf_norm) 
                                    last_inc = ic;
                                    interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                                    break
                                elseif dot(surf_norm,ray_dir) < acceptance_angle
                                    surface_collision_bool = true;   % collision with opaque surface voxel
                                    break
                                end
                        else
                            % Check refraction
                            nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                            if nn2 ~= nn1 % if refractive index changes
                                surface_collision_bool = false; % Medium refractive interface
                                last_inc = ic; % Last increment was ic
                                break
                            end
                        end
        
                        % Increment tau up to crossing XYZ(ia) boundary
                        tau_acc = tau_acc + (dist_frac_c)*default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                        if tau_acc > tau_ray 
                                                    % absorbed before crossing boundary
                            event = 2; % medium absorption
                            end_pos = XYZ;
                            rtrn = true;
                            break
                        end    
                    else % Hierarchy 2
                    %% did not increment XYZ(ib) or XYZ(ic)
        
                        % Increment tau up to crossing XYZ(ia) boundary
                        tau_acc = tau_acc + default_dist*VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)); 
                        if tau_acc > tau_ray 
                                                % absorbed before crossing boundary
                            event = 2; % medium absorption
                            end_pos = XYZ;
                            rtrn = true;
                            break
                        end
                    end
                end % Hierarchy 1
                %% Now Increment XYZ(ia)
                
                % Increment XYZ(ia)
                XYZ(ia) = XYZ(ia) + sgn(ia);      
                
                % Check exit or collision
                if XYZ(ia) > size_VS(ia) 
                                            % exited at maximum
                    if reflective_BCs(2,ia)     
                        sgn(ia) = -sgn(ia);         % flip direction
                        ray_dir(ia) = -ray_dir(ia);
                        XYZ(ia) = size_VS(ia);    % undo increment
                    else
                        end_pos = XYZ;            % note that X will outside of voxel_space
                        event = 1; % exit
                        rtrn = true;
                        break
                    end
                elseif XYZ(ia) < 1 
                                            % exited at minimum
                    if reflective_BCs(1,ia)     
                        sgn(ia) = -sgn(ia);         % flip direction
                        ray_dir(ia) = -ray_dir(ia);
                        XYZ(ia) = 1;              % undo increment
                    else
                        end_pos = XYZ;            % note that X will outside of voxel_space
                        event = 1; % exit
                        rtrn = true;
                        break
                    end
                elseif VS_opaq(XYZ(1),XYZ(2),XYZ(3))
                    surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                    if isempty(surf_norm) 
                        last_inc = ia;
                        interior_collision_bool = true; % edge case, crossed into non-surface opaque voxel
                        break
                    elseif dot(surf_norm,ray_dir) < acceptance_angle
                        surface_collision_bool = true;   % collision with opaque surface voxel
                        break
                    end
                else
                    % Check refraction
                    nn2 = VS_nn(XYZ(1),XYZ(2),XYZ(3));
                    if nn2 ~= nn1 % if refractive index changes
                        surface_collision_bool = false; % Medium refractive interface
                        last_inc = ia; % Last increment was ia
                        break
                    end
                end
                % Increment pb and pc
                pb = pb + Db;
                pc = pc + Dc;
            end % Traversal loop 1: If exited, the ray has collided with an opaque surface, reached a medium interface, or exited the domain
            if rtrn
                break
            end
        
            %%%%%%%%%%% EXITED TRAVERSAL %%%%%%%%%%%%%%%
            if surface_collision_bool
                %% Collision with opaque surface
                if rand <= VS_opaq_eps(XYZ(1),XYZ(2),XYZ(3))
                    %% absorbed by surface
                    end_pos = XYZ;
                    event = 3; % surface absorption
                    break
                else %% Currently only diffuse reflection
                    %% reflected by surface            
                    % Get surface norm
                    surf_norm = VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)};
                    %% Generate new ray    
                    % Generate emission direction (lambertian emitter)
                    sin_theta = sqrt(rand);                 % Modest Eq. 8.42
                    phi=2*pi*rand;                          % Modest Eq. 8.41
                    % Local direction vector using relation cos(asin(X)) = sqrt(1-X^2)
                    dir_local = [sin_theta.*cos(phi), sin_theta.*sin(phi), sqrt(1-sin_theta.^2)]; 
                    
                    % Transform to global coord (same logic as TransfoormCoord.m, but inlined for speed, and a slightly 
                    % different implementation by using the identity (AB)'=B'A' (i.e., (Rv')'=vR
                    % Assumes that surf_norm is already unit-length
                    if abs(surf_norm(3) - 1) < 10^(-6) % already [0 0 1]
                        ray_dir = dir_local;
            
                    elseif abs(surf_norm(3) + 1) < 10^(-6) % already [0 0 -1]
                        ray_dir = [-dir_local(1),dir_local(2),-dir_local(3)]; % rotate 180 degrees about y axis
            
                    else
                        nxy = sqrt(1-surf_norm(3)^2);   % (=== sin(theta) of normal vector)
                                                % Get new direction (matrix is equivalent to R' from TransformCoord.m)
                        ray_dir = dir_local*[ surf_norm(3)*surf_norm(1)/nxy, surf_norm(3)*surf_norm(2)/nxy,-nxy;...
                                             -surf_norm(2)/nxy,              surf_norm(1)/nxy,              0;...
                                              surf_norm(1),                  surf_norm(2),                  surf_norm(3)];
                    end
                    
                    % Generate new ray position
                    dist = 1/2*1/max(abs(surf_norm))-1e-10;   % Get distance to border of voxel from center of voxel.
                                                % New position in 1x1 plane normal to surface vector and centered on the
                                                % intersection of the surface normal and voxel boundary. 
                    ray_pos = (XYZ - 0.5)+ dist*surf_norm;
                end
            else % Not a surface collision
                %% Resolve Medium Interface or Edge Case of opaque non-surface
                % Initianlize bool
                specular_reflection_bool = false;
        
                % Figure out ray dir at collision (apply correction) since sign may have changed due to specular reflections on boundaries
                sgn_prev = sign(ray_dir);
                for ii = 1:3
                    if sgn_prev(ii) ~= sgn(ii)
                        ray_dir(ii) = -ray_dir(ii);
                    end
                end
        
                % Get position before increment (we use a new variable since we still need to know what XYZ is
                XYZ_prev = XYZ;
                XYZ_prev(last_inc) = XYZ_prev(last_inc) - sgn(last_inc);

                surf_norm = VS_surf_norms{XYZ_prev(1),XYZ_prev(2),XYZ_prev(3)};
                if interior_collision_bool % Just reemit from previous voxel
                    % This is the most problematic part of the code, in general
                    % If the domain is not nice, we can get infinite loops of reemisions (up to the maximum loop_count)
                    % These stuck rays take up a lot of computation time and don't benefit the solution
                    if dot(surf_norm,ray_dir) < acceptance_angle
                        end_pos = XYZ_prev;
                        event = 3;
                        break
                    end
                    % Generate new ray direction from previous XYZ_prev
                    dist = 1/2*1/max(abs(surf_norm));
                    XYZ_next = floor(XYZ_prev-0.5 + (dist+1e-10)*surf_norm)+1;
                    if any(XYZ_next == 0) || any(XYZ_next > size_VS)
                        event = 1;
                        end_pos = XYZ_next;
                        break;
                    end
                    if VS_opaq(XYZ_next(1),XYZ_next(2),XYZ_next(3)) && isempty(VS_surf_norms(XYZ_next(1),XYZ_next(2),XYZ_next(3)))
                        disp("Ray is trapped!")
                        event = 3;
                        end_pos = XYZ_prev;
                        break;
                    end
                    ray_pos = (XYZ_prev - 0.5)+ (dist-1e-10)*surf_norm;
                    
                    sin_theta = sqrt(rand);                 % Modest Eq. 8.42
                    phi=2*pi*rand;                          % Modest Eq. 8.41
                    dir_local = [sin_theta.*cos(phi), sin_theta.*sin(phi), sqrt(1-sin_theta.^2)]; 
                    if abs(surf_norm(3) - 1) < 10^(-6) % already [0 0 1]
                        ray_dir = dir_local;
                    elseif abs(surf_norm(3) + 1) < 10^(-6) % already [0 0 -1]
                        ray_dir = [-dir_local(1),dir_local(2),-dir_local(3)]; % rotate 180 degrees about y axis
            
                    else
                        nxy = sqrt(1-surf_norm(3)^2);   % (=== sin(theta) of normal vector)
                                                % Get new direction (matrix is equivalent to R' from TransformCoord.m)
                        ray_dir = dir_local*[ surf_norm(3)*surf_norm(1)/nxy, surf_norm(3)*surf_norm(2)/nxy,-nxy;...
                                             -surf_norm(2)/nxy,              surf_norm(1)/nxy,              0;...
                                              surf_norm(1),                  surf_norm(2),                  surf_norm(3)];
                    end
                else % Medium refraction
                
                    % Determine exact position of interface using values of p2 and p3 and last_inc
                    if sgn(last_inc) == 1
                        ray_pos(last_inc) = XYZ_prev(last_inc);
                    else % Sign is -1
                        ray_pos(last_inc) = XYZ_prev(last_inc) - 1;
                    end
            
                    if last_inc == ia % If driving axis incremented last pb and pc tell us the distance from voxel boundary
                        if sgn(ib) == 1
                            ray_pos(ib) = XYZ_prev(ib) + pb;
                        else % Sign is -1
                            ray_pos(ib) = XYZ_prev(ib) - pb - 1;
                        end
                        
                        if sgn(ic) == 1
                            ray_pos(ic) = XYZ_prev(ic) + pc;
                        else
                            ray_pos(ic) = XYZ_prev(ic) - pc - 1;
                        end
                    elseif last_inc == ib % pb > 0 -> determine coordinates with dist_frac_b and pc < c
                        if sgn(ia) == 1
                            ray_pos(ia) = XYZ_prev(ia) - dist_frac_b;
                        else
                            ray_pos(ia) = XYZ_prev(ia) + dist_frac_b - 1;
                        end
                        
                        if sgn(ic) == 1
                            ray_pos(ic) = XYZ_prev(ic) + pc - dist_frac_b*Dc; % pc gives position when crossing a-axis and then move backward by dist_frac_b*Dc to get position crossing b axis
                        else
                            ray_pos(ic) = XYZ_prev(ic) - pc + dist_frac_b*Dc - 1;
                        end
                    else % last_inc == ic -> determine coordinates with dist_frac_c and pb < 0
                        if sgn(ia) == 1
                            ray_pos(ia) = XYZ_prev(ia) - dist_frac_c;
                        else
                            ray_pos(ia) = XYZ_prev(ia) + dist_frac_c - 1;
                        end
            
                        if sgn(ib) == 1
                            ray_pos(ib) = XYZ_prev(ib) + pb - dist_frac_c*Db;
                        else
                            ray_pos(ib) = XYZ_prev(ib) - pb + dist_frac_c*Db - 1;
                        end
                    end
                    % With ray position determined, we can determine the new ray direction after the medium collision
                    
                    if isempty(surf_norm) % We will just use the voxel boundary as the plane (this should never occur, though)
                        surf_norm = zeros(1,3);
                        surf_norm(last_inc) = sgn(last_inc); % Points in direction such that dot(surf_norm,ray_dir)>0
                    end
                    
                    % cos(theta) = dot(a,b)/(|a|x|b|): https://chortle.ccsu.edu/vectorlessons/vch07/vch07_3.html
                    cos1 = dot(ray_dir,surf_norm); % Both vectors are already magnitude 1 and take negative since we want acute angle < pi/2
                    adjust_pos = false;
                    if cos1 < 0 % edge case where cos1 < 0 can occur due to discretization of voxel space
                        adjust_pos = true; % Boolean used later
                    else
                        theta1 = acos(cos1);
                        mu = nn1/nn2; % Ratio of refractive indices
                        
                        if mu > 1 % Going from larger refractive index to smaller
                            % Check for total internal reflection
                            theta_c = asin(1/mu);
                            if theta1 > theta_c
                                % Specular reflection (resolved a few lines later)
                                specular_reflection_bool = true;
                            end
                        end
            
                        % Maybe reflection, maybe refraction
                        % Need cosine of each angle
                        % Already have cos1
                        if ~specular_reflection_bool
                            cos2 = cos(asin(sin(theta1)*mu)); % From Snell's law n1sin(theta1)=n2sin(theta2)
                            
                            rho = 1/2*(((nn1*cos2-nn2*cos1)/(nn1*cos2+nn2*cos1)).^2+((nn1*cos1-nn2*cos2)/(nn1*cos1+nn2*cos2)).^2); % Modest Eq: 2.96
                
                            if rand() < rho % Ray is reflected, since if rho == 1, should always be reflected
                                % Reflect specularly (resolved a few lines later)
                                specular_reflection_bool = true;
                            else
                                % Transmit and refract 
                                ray_dir = mu*ray_dir + surf_norm*sqrt(1-mu^2*(1-cos1^2)) - mu*surf_norm*cos1;
                                sgn = sign(ray_dir);
                                % https://opg.optica.org/josaa/abstract.cfm?uri=josaa-29-7-1358 Equation (3)
                                % in their notation:
                                %   g = surf_norm
                                %   s' = ray_dir_new
                                %   s = ray_dir,
                                %   (gs) = dot(g,s) = -dot(ray_dir,surf_norm) = cos1
                            end
                        end
                        
                        if specular_reflection_bool % Resolve specular reflection
                            ray_dir = ray_dir - 2*cos1*surf_norm; % Update ray direction by flipping component along surf_norm
                            sgn_new = sign(ray_dir);
                            % We need to check for edge case where there is a specular reflection yet the ray still ends up in
                            % the refracting voxel. This is possible since surf_norm is often not normal to the voxel face
                            % e.g., for curved surfaces, inclined planes. This edge case may even be quite common    
                            if sign(ray_dir(last_inc)) == sgn(last_inc) % We will end up in medium 2 voxel despite specular reflection 
                               adjust_pos = true; 
                            end
                            sgn = sgn_new;
                        end
                    end
                    if adjust_pos % Edge case has occurred
                        % We want to avoid crossing the boundary
                        dist_to_clear = (0.5 - (sgn.*(rem(ray_pos,1)-0.5))).*vx_scale; % Get distance to next boundary in direction of ray travel (note this is the opposite of xs calculated in Premable pt 2)
                        time_to_clear = dist_to_clear./abs(ray_dir); % "time" because direction components are essentially velocities
                        i_rem = [1,2,3];
                        i_rem(last_inc) = [];
                        t_min = inf;
                        for ii = i_rem 
                            if time_to_clear(ii) < t_min
                                i_min = ii;
                                t_min = time_to_clear(ii);
                            end
                        end
                        ray_pos(last_inc) = ray_pos(last_inc) - 0.1*sgn(last_inc); % Move ray a little bit backward from boundary
                        ray_pos(i_min) = ray_pos(i_min)+dist_to_clear(i_min)/vx_scale(i_min)*(1-1e-10); % Move ray so that the next boundary it crosses will be i_min
                    end                
                end % end: interior collision bool
                % IT TURNS OUT THIS NEXT LINE IS EXTREMELY IMPORTANT
                ray_pos = ray_pos + sgn*eps.*XYZ; % Give small bump to get off boundary (scale bump with position for floating point reasons (floor(1-eps) = 0, but floor(10-eps) = 10)
            end   % end surface collision bool  
        %% Go back to Preamble part 2
        end % Preamble loop
        if loop_count == max_loops && event ~= 1 % Resolve stuck rays (as best we can!)
            if VS_opaq(XYZ(1),XYZ(2),XYZ(3)) % Check current voxel
                if VS_opaq_eps(XYZ(1),XYZ(2),XYZ(3))>0
                    event = 3;
                    end_pos = XYZ;
                    rtrn = true;
                else % If it's a pure reflective voxel -> move it along the surface normal into the next voxel
                    XYZ = ceil(XYZ-0.5 + VS_surf_norms{XYZ(1),XYZ(2),XYZ(3)});
                end
            end
            if VS_PM_kappa(XYZ(1),XYZ(2),XYZ(3)) > 0 && ~rtrn % Check participating media
                event = 2;
                end_pos = XYZ;
                rtrn = true;
            end
            if ~rtrn % Pick a random absorption location (might be expensive, but hopefully this is not called frequently)
                locs = find(VS_opaq_eps>0 | VS_PM_kappa>0);
                [X,Y,Z] = ind2sub(size_VS,locs(randi(length(locs))));
                end_pos = [X,Y,Z];
                if VS_opaq_eps(X,Y,Z)>0
                    event = 3;
                else
                    event = 2;
                end
            end
        elseif event == 2 && all(end_pos == XYZ_0) % Self absorption in participating media voxel
            event = 4;
        end
        rays_end_pos(nn_ray,:) = end_pos;
        rays_events(nn_ray) = event;
        
    end % End parfor
end % End function