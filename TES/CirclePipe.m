classdef CirclePipe < PipeSystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius_vx;
    end
    
    methods
        function obj = CirclePipe(x_loc,y_loc,radius_vx,roughness,vx_scale,size_VS,fluid_property_table)
            %PIPESYSTEM Construct an instance of this class
            %   Detailed explanation goes here
           
            obj.x_loc = x_loc;
            obj.y_loc = y_loc;
            obj.radius_vx = radius_vx;
            obj.roughness = roughness;
            obj.vx_scale = vx_scale;
            obj.size_VS = size_VS;
            obj.D_hydraulic = radius_vx*2*vx_scale(1);
            obj.A_cs = pi/4*obj.D_hydraulic^2;
            [VS_pipe_2D,VS_boundaries_2D,inds_pipe_boundary,inds_pipe_interior,boundary_perimiter_span] =  getPipeGeometry(obj);
            obj.VS_pipe_2D = VS_pipe_2D;
            obj.VS_boundaries_2D = VS_boundaries_2D;
            obj.inds_pipe_boundary = inds_pipe_boundary;
            obj.inds_pipe_interior = inds_pipe_interior;
            obj.boundary_perimiter_span = boundary_perimiter_span;
            obj.fluid_property_table = fluid_property_table;
            obj.S_wall_z = zeros(size_VS(3),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getPipeGeometry
        function [VS_pipe,VS_boundaries,inds_pipe_boundary,inds_pipe_interior,boundary_perimiter_span] =  getPipeGeometry(obj,plot_bool)
            if nargin == 1
                plot_bool = false;
            end
            ri = obj.radius_vx; x_center = obj.x_loc; y_center = obj.y_loc;
            Nx = obj.size_VS(1);Ny = obj.size_VS(2);
            
            ii_center = floor(x_center)+1;
            jj_center = floor(y_center)+1;
            count = 0;
            inds_pipe = [];
            VS_pipe = false(Nx,Ny);
            for dx_vx = round(-(ri+1)):round(ri+1)
                ii = ii_center + dx_vx;
                dx_cur = ii-0.5-x_center;
                for dy_vx = round(-(ri+1)):round(ri+1)
                    jj = jj_center + dy_vx;
                    dy_cur = jj-0.5-y_center;
                    r_cur = (dx_cur^2 + dy_cur^2)^0.5;
                    if r_cur <= ri
                        count = count + 1;
                        inds_pipe(count,:) = [ii,jj];
                        VS_pipe(ii,jj) = true;
                    end
                end
            end
            N_pipe_cells = size(inds_pipe,1);
            VS_boundaries = 6-countConnectivity(VS_pipe);
            VS_boundaries(~VS_pipe) = 0;
            N_boundary = sum(VS_boundaries(:)>0);
            
            if N_pipe_cells == 1
                inds_pipe_boundary = inds_pipe;
                inds_pipe_interior = [];
                boundary_perimiter_span = 2*ri*pi;
            else
                thetas_cell_boundary = zeros(N_boundary,1);
                inds_pipe_boundary = zeros(N_boundary,2);
                inds_pipe_interior = zeros(N_pipe_cells-N_boundary,2);
                count_boundary = 0;
                count_interior = 0;
                for nn = 1:N_pipe_cells
                    ii = inds_pipe(nn,1);
                    jj = inds_pipe(nn,2);
                    if VS_boundaries(ii,jj)>0 % Then it is a boundary cell
                        count_boundary = count_boundary+1;
                        inds_pipe_boundary(count_boundary,:) = [ii,jj];
                        thetas_cell_boundary(count_boundary) = atan2(jj-0.5-y_center,ii-0.5-x_center);
                    else
                        count_interior = count_interior+1;
                        inds_pipe_interior(count_interior,:) = [ii,jj];
                    end
                end
                [thetas_cell_boundary,sort_inds] = sort(thetas_cell_boundary);
                inds_pipe_boundary = inds_pipe_boundary(sort_inds,:);
                thetas_cell_boundary(end+1) = thetas_cell_boundary(1)+2*pi; % Make periodic
                diff_thetas = diff(thetas_cell_boundary); % Take difference
                thetas_cell_boundary(end) = [];
                thetas_cell_boundary = [thetas_cell_boundary(end)-2*pi;thetas_cell_boundary];
                diff_thetas = [diff(thetas_cell_boundary),diff_thetas];
                thetas_cell_boundary(1) = [];
                
                thetas_cell_edges = zeros(N_boundary,2);
                for nn = 1:N_boundary
                    theta_cell = thetas_cell_boundary(nn);
                    ii = inds_pipe_boundary(nn,1);
                    jj = inds_pipe_boundary(nn,2);
                
                    theta_b = theta_cell;
                    theta_a = theta_b - 0.99*diff_thetas(nn,1);
                    err = 1e6;
                    while err>1e-6
                        theta_mid = (theta_a + theta_b)/2;
                        p_mid = [x_center,y_center] + [cos(theta_mid),sin(theta_mid)]*ri*0.9999999;
                        ii_mid = floor(p_mid(1))+1;
                        jj_mid = floor(p_mid(2))+1;
                        [~,max_vec_ind] = max(abs([cos(theta_mid),sin(theta_mid)]));
                        if max_vec_ind == 1 % x-component larger -> consider y boundaries
                            vx_dif = abs(jj_mid-jj);
                        else % consider x- boundaries
                            vx_dif = abs(ii_mid-ii); 
                        end
                        if vx_dif > 0
                            theta_a = theta_mid;
                        else
                            theta_b = theta_mid;
                        end
                        if vx_dif > 1
                            err = 1e6;
                        else
                            err = min(abs(rem(p_mid,1)));
                        end
                    end
                    thetas_cell_edges(nn,1) = theta_mid;
                    theta_a = theta_cell;
                    theta_b = theta_a + 0.99*diff_thetas(nn,2);
                    err = 1e6;
                    while err>1e-6
                        theta_mid = (theta_a + theta_b)/2;
                        p_mid = [x_center,y_center] + [cos(theta_mid),sin(theta_mid)]*ri*0.9999999; % Adjustment by factor of 0.9999999 prevents edge-cases on boundaries
                        ii_mid = floor(p_mid(1))+1;
                        jj_mid = floor(p_mid(2))+1;
                        [~,max_vec_ind] = max(abs([cos(theta_mid),sin(theta_mid)]));
                        if max_vec_ind == 1 % x-component larger -> consider y boundaries
                            vx_dif = abs(jj_mid-jj);
                        else % consider x- boundaries
                            vx_dif = abs(ii_mid-ii);     
                        end
                        if vx_dif > 0
                            theta_b = theta_mid;
                        else
                            theta_a = theta_mid;
                        end
                        if vx_dif > 1
                            err = 1e6;
                        else
                            err = min(abs(rem(p_mid,1)));
                        end
                    end
                    thetas_cell_edges(nn,2) = theta_mid;
                end
                thetas_cell_edges_adj = zeros(N_boundary,2);
                for nn = 1:N_boundary
                    nn_next = mod(nn,N_boundary)+1; % Circular indexing
                    if abs(thetas_cell_edges(nn,2)-thetas_cell_edges(nn_next,1))>pi % They are separated by 2pi
                        theta_l = min(thetas_cell_edges(nn,2),thetas_cell_edges(nn_next,1));
                        theta_u = max(thetas_cell_edges(nn,2),thetas_cell_edges(nn_next,1));
                        theta_shift = theta_l+2*pi;
                        theta_mid = (theta_shift+theta_u)/2;
                        if thetas_cell_edges(nn,2)<thetas_cell_edges(nn_next,1)
                            thetas_cell_edges_adj(nn,2) = theta_mid-2*pi;
                            thetas_cell_edges_adj(nn_next,1) = theta_mid;
                        else
                            thetas_cell_edges_adj(nn_next,1) = theta_mid-2*pi;
                            thetas_cell_edges_adj(nn,2) = theta_mid;
                        end
                    else
                        theta_mid = (thetas_cell_edges(nn,2) + thetas_cell_edges(nn_next,1))/2;
                        thetas_cell_edges_adj(nn,2) = theta_mid;
                        thetas_cell_edges_adj(nn_next,1) = theta_mid;
                    end
                end
                boundary_perimiter_span = (thetas_cell_edges_adj(:,2)-thetas_cell_edges_adj(:,1))*ri;
                if abs(sum(boundary_perimiter_span)-2*pi*ri) > 1e-6
                    error("Perimiter does not add up!")
                end
            end 
            
            if plot_bool
                figure;
                hold on;
                grid on;
                axis equal;
                for nn = 1:N_pipe_cells
                    set(gca,'Color',[0.8,0.8,0.8])
                    x_fill = inds_pipe(nn,1);
                    y_fill = inds_pipe(nn,2);
                    if ismember([x_fill,y_fill],inds_pipe_boundary,'rows')
                        color = [0.9,0.9,0.9];
                    else
                        color = [1,1,1];
            
                    end
                    x_fill = [x_fill,x_fill,x_fill-1,x_fill-1];
                    y_fill = [y_fill,y_fill-1,y_fill-1,y_fill];
                    fill(x_fill,y_fill,color,'FaceAlpha',1,'LineWidth',2)  
                end
                if N_boundary ~= 1
                    x_cell_edges = x_center + cos(thetas_cell_edges(:))*ri;
                    y_cell_edges = y_center + sin(thetas_cell_edges(:))*ri;
                    
                    x_cell_edges_adj = x_center + cos(thetas_cell_edges_adj(:))*ri;
                    y_cell_edges_adj = y_center + sin(thetas_cell_edges_adj(:))*ri;
                    
                    scatter(x_cell_edges_adj,y_cell_edges_adj,'o','MarkerEdgeColor','k','SizeData',20)
                    if not(all(abs(x_cell_edges-x_cell_edges_adj)<1e-3) && all(abs(y_cell_edges-y_cell_edges_adj)<1e-3))
                        scatter(x_cell_edges_adj,y_cell_edges_adj,'x','MarkerEdgeColor','k','MarkerFaceColor','k')
                    end
                end
                all_thetas = linspace(-pi,pi,1000);
                all_x = x_center + cos(all_thetas)*ri;
                all_y = y_center + sin(all_thetas)*ri;
                plot(all_x,all_y,'Color',	"#0072BD",'LineWidth',1);
                xlim([x_center-round(ri)-1,x_center+round(ri)+1])
                ylim([y_center-round(ri)-1,y_center+round(ri)+1])
                xtick_bounds = xlim();
                xtick_bounds(1) = ceil(xtick_bounds(1));
                xtick_bounds(2) = floor(xtick_bounds(2));
                ytick_bounds = ylim();
                ytick_bounds(1) = ceil(ytick_bounds(1));
                ytick_bounds(2) = floor(ytick_bounds(2));
                xticks(xtick_bounds(1):1:xtick_bounds(2));
                yticks(ytick_bounds(1):1:ytick_bounds(2));
                set(gca,'GridColor',[0 0 0]);
                set(gca,'GridLineWidth',2)
               
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getPipeEquation
        
        function [M_pipe, RHS_pipe] = getPipeEquation(obj,mesh,k,h,T_b)
            % GETPIPEEQUATION Define a convective boundary condition at index ii which gives the heat transfer equivalent to a cylindrical
            % pipe (its axis orthogonal to the x-axis) which has heat transfer coefficient h through a solid with heat conductivity k
            % It is assumed that flow is symmetric along the perimiter of the cylinder and heat conductivity is the same in both
            % boundary cells. Note that this assumes the uniform cartesian grid!
            % Credit: Eftekhari, A.A. et al. (2015). FVTool: a finite volume toolbox for Matlab. Zenodo. http://doi.org/10.5281/zenodo.32745
            % Parameters:
            %   mesh: Domain mesh variable
            %   k [W/(K-m)]: (scalar) Heat conductivity of the wall material
            %   h [W/(K-m^2))]: (scalar) Heat transfer coefficient fluid at each z-location along the pipe axis
            %   T_b [K]: (1D vector) Fluid temperature of fluid along pipe
            %   S [K/s]: Heat generation in boundary voxels along axis
            S = obj.S_wall_z;
            Nxyz = mesh.dims;
            if length(Nxyz) == 2
                Nxyz(3) = 1;
            end
            Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
            x_center = obj.x_loc; y_center = obj.y_loc; ri = obj.radius_vx; % I think this is a little more readable

            ii_center = floor(x_center)+1;
            jj_center = floor(y_center)+1;
            
            
            G = reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
            dx = mesh.cellsize.x(ii_center+1);
            dy = mesh.cellsize.y(jj_center+1);
            
            if dx ~= dy
                disp("Warning, nonuninform mesh (getPipeEquation)")
            end
            dxy = dx;
            
            % define the vectors to be used for the creation of the sparse matrix
            N_boundary = size(obj.inds_pipe_boundary,1);
            N_interior = size(obj.inds_pipe_interior,1);
            nb = (N_boundary*5+N_interior)*Nz;
            iii = zeros(nb,1);
            jjj = zeros(nb,1);
            s = zeros(nb,1);
            
            % define the RHS column vector
            RHS_pipe = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
            
            % Define matrix of coefficients and RHS for ghost cell for solving system
            q = 0;
            % interior cells -> have no information set temperature to 0;
            for nn = 1:N_interior
                
                ii = obj.inds_pipe_interior(nn,1)+1;
                jj = obj.inds_pipe_interior(nn,2)+1;
                q=q(end)+(1:Nz);
                iii(q) = G(ii,jj,2:end-1);  jjj(q) = G(ii,jj,2:end-1); s(q) = 1; % at pipe
            end
            for nn = 1:N_boundary
                ii = obj.inds_pipe_boundary(nn,1)+1;
                jj = obj.inds_pipe_boundary(nn,2)+1;
                perimiter = obj.boundary_perimiter_span(nn);
                N_faces = obj.VS_boundaries_2D(ii-1,jj-1); % Number of faces the cell has which comprise the interface
                h_star = h*perimiter/N_faces; % modified heat transfer coefficient based on the ratio of the "actual" surface area of the cylidner to the surface area of the discretized interface.
                
                iii_pipe = squeeze(G(ii,jj,2:end-1)); iii_left = squeeze(G(ii-1,jj,2:end-1)); iii_right = squeeze(G(ii+1,jj,2:end-1)); iii_bot = squeeze(G(ii,jj-1,2:end-1)); iii_top = squeeze(G(ii,jj+1,2:end-1));
                % At cell location
                K1 = k/dxy;
                q = q(end)+(1:Nz);
                iii(q) = iii_pipe;  jjj(q) = iii_pipe;  s(q) = K1; % At cell
                
                % check left boundary
                K3 = 0;
                if ~obj.VS_pipe_2D(ii-2,jj-1)
                    ro = ((ii-2-0.5-x_center).^2+(jj-1-0.5-y_center).^2)^0.5;
                    K_const = 1-log(ro/ri)./(log(ro/ri)+k./h/dxy/ri); % Note regular h in denominator, not h_star!
                    % These values obtained from reformulating a boundary condition of the third kind into a*dT/dx+b*T=c
                    K2 = -k/dxy + h_star.*K_const; 
                    K3 = K3 + (h_star.*T_b.*K_const+S.*h_star/4./k*(ro^2-ri^2)*dxy^2.*K_const+S.*h_star/2./h.*(1-K_const)); % K3 is an average of K3 obtained for each interface
                else
                    K2 = 0;
                end
                q=q(end)+(1:Nz);
                iii(q) = iii_pipe;  jjj(q) = iii_left;  s(q) = K2(:)/N_faces; % left of pipe.
                % check right boundary
                if ~obj.VS_pipe_2D(ii,jj-1)
                    ro = ((ii-0.5-x_center).^2+(jj-1-0.5-y_center).^2)^0.5;
                    K_const = 1-log(ro/ri)./(log(ro/ri)+k./h/dxy/ri); % Note regular h in denominator, not h_star!
                    K2 = -k/dxy + h_star.*K_const; 
                    K3 = K3 + (h_star.*T_b.*K_const+S.*h_star/4./k*(ro^2-ri^2)*dxy^2.*K_const+S.*h_star/2./h.*(1-K_const)); % K3 is an average of K3 obtained for each interface
                else
                    K2 = 0;
                end
                q=q(end)+(1:Nz);
                iii(q) = iii_pipe; jjj(q) = iii_right; s(q) = K2/N_faces; % right of pipe
                % check bottom boundary
                if ~obj.VS_pipe_2D(ii-1,jj-2)
                    ro = ((ii-1-0.5-x_center).^2+(jj-2-0.5-y_center).^2)^0.5;
                    K_const = 1-log(ro/ri)./(log(ro/ri)+k./h/dxy/ri); % Note regular h in denominator, not h_star!
                    K2 = -k/dxy + h_star.*K_const; 
                    K3 = K3 + (h_star.*T_b.*K_const+S.*h_star/4./k*(ro^2-ri^2)*dxy^2.*K_const+S.*h_star/2./h.*(1-K_const)); % K3 is an average of K3 obtained for each interface
                else
                    K2 = 0;
                end
                q = q(end)+(1:Nz);
                iii(q) = iii_pipe;  jjj(q) = iii_bot;  s(q) = K2/N_faces; % bottom of pipe
            
                % check top boundary
                if ~obj.VS_pipe_2D(ii-1,jj)
                    ro = ((ii-1-0.5-x_center).^2+(jj-0.5-y_center).^2)^0.5;
                    K_const = 1-log(ro/ri)./(log(ro/ri)+k./h/dxy/ri); % Note regular h in denominator, not h_star!
                    K2 = -k/dxy + h_star.*K_const; 
                    K3 = K3 + (h_star.*T_b.*K_const+S.*h_star/4./k*(ro^2-ri^2)*dxy^2.*K_const+S.*h_star/2./h.*(1-K_const)); % K3 is an average of K3 obtained for each interface
                else
                    K2 = 0;
                end
                q=q(end)+(1:Nz);
                iii(q) = iii_pipe; jjj(q) = iii_top; s(q) = K2/N_faces; % top of pipe
                RHS_pipe(iii_pipe) = K3/N_faces;
            end
            % Build the sparse matrix of the boundary conditions
            M_pipe = sparse(iii, jjj, s, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));  
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getWallTemperature
        function [T_wall_z] = getWallTemperature(obj,VS_T,thermal_conductivity,h_mean_z,T_fluid_z)
            TES_boundary_inds = find((6 - countConnectivity(~obj.VS_pipe_2D)).*(~obj.VS_pipe_2D));
            [ii,jj] = ind2sub(size(obj.VS_pipe_2D),TES_boundary_inds);
            N_TES_boundary = length(TES_boundary_inds);
            T_wall_all = zeros(obj.size_VS(3),N_TES_boundary);
            dxy = obj.vx_scale(1);
            ri = obj.radius_vx;
            for nn = 1:N_TES_boundary
                T_cur = squeeze(VS_T(ii(nn),jj(nn),:)); % Column vector
                ro = ((ii(nn)-0.5-obj.x_loc)^2 + ...
                            (jj(nn)-0.5-obj.y_loc)^2).^0.5;
                if isscalar(thermal_conductivity)
                    k = thermal_conductivity;
                else
                    k = squeeze(thermal_conductivity(ii,jj,:));
                end
                %T_wall_all(:,nn) = T_cur - (T_cur-T_fluid_z)./(k./h_mean_z./(ri*dxy)+log(ro/ri))*log(ro/ri);
                K_const = 1-log(ro/ri)./(log(ro/ri)+k./h_mean_z/dxy/ri);
                T_wall_all(:,nn) = T_cur +obj.S_wall_z/4./k*(ro.^2-ri^2)*dxy^2.*K_const + obj.S_wall_z/2./h_mean_z*ri*dxy.*(1-K_const)+ (T_fluid_z-T_cur).*(1-K_const);
            end
            T_wall_z = mean(T_wall_all,2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% get_f_Darcy
        function f_Darcy = get_f_Darcy(obj,Re,roughness)
            %% Compute darcy friction factor
            % Re: Reynolds number
            % using the weighting method given by Avci (2019), but using the Haaland correlation for the turbulent friction factor
            rel_roughness = roughness/(obj.D_hydraulic); % Relative roughness = epsilon/D
            
            if Re < 4000
                if Re<=2300
                    Re_lam = Re;
                else
                    Re_lam = 2300;
                end
                f_lam = 64/Re;
            end
            if Re > 2300
                if Re >= 4000
                    Re_turb = Re;
                else
                    Re_turb = 4000;
                end
                f_t = (-1.8*log10((rel_roughness/3.7).^1.11+6.9./Re)).^(-2); % Haaland correlation for turbulent region
            end
            if Re <= 2300
                f_Darcy = f_lam;
            elseif Re >= 4000
                f_Darcy = f_t;
            else
                alpha = (Re-2300)/(4000-2300); % Simple linear interpolation of f_lam and f_t;
                f_Darcy = (1-alpha)*f_lam + alpha*f_t;
            end
        end
    end
    methods (Static)
        function Nu_x = getNu_x_lam(x_star)
            % Gz = Graetz number
            % Lienhard, Eq 7.28, or Bhatti and Shah, Eq 3.48 (note misprint, condition is x<= 0.001, not x<=0.01
            if x_star <= 0.001 % 
                Nu_x = 1.077*x_star.^(-1/3)-0.7;
            else
                Nu_x = 3.657+0.2362.*x_star.^(-0.488).*exp(-57.2.*x_star);
            end
        end
        function Nu_m = getNu_m_lam(x_star_L,x_star_U)
            % x_star = 1/Gz
            % x_star_L < x_star_U
            if x_star_L == 0
                Nu_m = 3.657./tanh(2.264*x_star_U.^(1/3)+1.7*x_star_U.^(-2/3))+0.0499./x_star_U*tanh(x_star_U); % Lienhard Eq 7.29, replace Gz with 1/x_star_U
            elseif x_star_U <= 0.001
                Nu_m = integral(@(x) 1.077*x.^(-1/3)-0.7,x_star_L,x_star_R);
            elseif x_star_L >= 0.001
                Nu_m = integral(@(x) 3.657+0.2362.*x.^(-0.488).*exp(-57.2.*x),1/1000,x_star_U);
            else % x_star_L < 1/1000 && x_star_R > 1/1000
                Nu_m = integral(@(x) 1.077*x.^(-1/3)-0.7,x_star_L,1/1000) + integral(@(x) 3.657+0.2362.*x.^(-0.488).*exp(-57.2.*x),1/1000,x_star_U);
            end
            Nu_m = Nu_m/(x_star_U-x_star_L);
        end
    end
end

