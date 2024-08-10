classdef SquarePipe < PipeSystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xy_width_vx; % (scalar) voxel units
    end
    
    methods
        function obj = SquarePipe(x_loc,y_loc,xy_width_vx,roughness,vx_scale,size_VS,fluid_property_table)
            %PIPESYSTEM Construct an instance of this class
            %   Detailed explanation goes here
           
            obj.x_loc = x_loc;
            obj.y_loc = y_loc;
            obj.xy_width_vx = xy_width_vx;
            obj.D_hydraulic = xy_width_vx*vx_scale(1);
            obj.A_cs = xy_width_vx^2*vx_scale(1)^2;
            obj.roughness = roughness;
            obj.vx_scale = vx_scale;
            obj.size_VS = size_VS;
            [VS_pipe_2D,VS_boundaries_2D,inds_pipe_boundary,inds_pipe_interior,boundary_perimiter_span] =  getPipeGeometry(obj);
            obj.VS_pipe_2D = VS_pipe_2D;
            obj.VS_boundaries_2D = VS_boundaries_2D;
            obj.inds_pipe_boundary = inds_pipe_boundary;
            obj.inds_pipe_interior = inds_pipe_interior;
            obj.boundary_perimiter_span = boundary_perimiter_span;
            obj.fluid_property_table = fluid_property_table;
            obj.S_wall_z = zeros(size_VS(3),1);
        end
        
        function [VS_pipe,VS_boundaries,inds_pipe_boundary,inds_pipe_interior,boundary_perimiter_span] =  getPipeGeometry(obj,plot_bool)
            if nargin == 1
                plot_bool = false;
            end
            xyw = obj.xy_width_vx; x_center = obj.x_loc; y_center = obj.y_loc;
            xl = x_center-xyw/2; xu = x_center+xyw/2; yl = y_center-xyw/2; yu = y_center+xyw/2;
            x_ind_range = [floor(xl+0.5)+1,floor(xu-0.5)+1];
            y_ind_range = [floor(yl+0.5)+1,floor(yu-0.5)+1];

            Nx = obj.size_VS(1);Ny = obj.size_VS(2);
            
            VS_pipe = false(Nx,Ny);
            VS_pipe(x_ind_range(1):x_ind_range(2),y_ind_range(1):y_ind_range(2)) = true;
            [X,Y] = meshgrid(x_ind_range(1):x_ind_range(2),y_ind_range(1):y_ind_range(2));
            inds_pipe = [X(:),Y(:)];

            N_pipe_cells = size(inds_pipe,1);
            VS_boundaries = 6-countConnectivity(VS_pipe);
            VS_boundaries(~VS_pipe) = 0;
            N_boundary = sum(VS_boundaries(:)>0);
           
            if N_pipe_cells == 1
                inds_pipe_boundary = inds_pipe;
                inds_pipe_interior = [];
                boundary_perimiter_span = 2*xyw+2*xyw;
            else
                inds_pipe_boundary = inds_pipe(inds_pipe(:,1) == x_ind_range(1) | inds_pipe(:,1) == x_ind_range(2) |...
                                          inds_pipe(:,2) == y_ind_range(1) | inds_pipe(:,2) == y_ind_range(2),:);
                inds_pipe_interior = inds_pipe(~(inds_pipe(:,1) == x_ind_range(1) | inds_pipe(:,1) == x_ind_range(2) |...
                                          inds_pipe(:,2) == y_ind_range(1) | inds_pipe(:,2) == y_ind_range(2)),:);
                boundary_perimiter_span = zeros(N_boundary,1);
                for nn = 1:N_boundary
                    cur_loc = inds_pipe_boundary(nn,:);
                    if VS_boundaries(cur_loc(1),cur_loc(2)) == 1
                        % It's an edge
                        boundary_perimiter_span(nn) = 1; 
                       
                    else % It's a corner
                        cur_corner_dist = abs(cur_loc-0.5-[x_center,y_center])+0.5;
                        cur_diff = cur_corner_dist - [xyw/2,xyw/2]; % Accounting for pipe which is not aligned with grid
                        boundary_perimiter_span(nn) = 2-sum(cur_diff);
                    end
                end
                if abs(sum(boundary_perimiter_span)-2*xyw-2*xyw) > 1e-6
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
                x_cell_edges = [x_center-xyw/2,x_center+xyw/2,x_center+xyw/2,x_center-xyw/2,x_center-xyw/2];
                y_cell_edges = [y_center-xyw/2,y_center-xyw/2,y_center+xyw/2,y_center+xyw/2,y_center-xyw/2];

                plot(x_cell_edges,y_cell_edges,'Color',	"#0072BD",'LineWidth',1);
                xlim([floor(x_center-xyw/2-2),ceil(x_center+xyw/2+2)])
                ylim([floor(y_center-xyw/2-2),ceil(y_center+xyw/2+2)])
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
        
        function [M_pipe, RHS_pipe] = getPipeEquation(obj,mesh,k,h,T_b)
            % GETPIPEEQUATION Define a convective boundary condition at index ii which gives the heat transfer equivalent to a cylindrical
            % pipe (its axis orthogonal to the x-axis) which has heat transfer coefficient h through a solid with heat conductivity k
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
            xyw = obj.xy_width_vx; x_center = obj.x_loc; y_center = obj.y_loc; % I think this is a little more readable
            xl = x_center-xyw/2; xu = x_center+xyw/2; yl = y_center-xyw/2; yu = y_center+xyw/2;
            offset_xl = xl-round(xl); offset_xu = round(xu)-xu; offset_yl = yl-round(yl); offset_yu = round(yu)-yu;
            

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
            dw_xl = dxy*(1/2+offset_xl); dw_xu = dxy*(1/2+offset_xu); dw_yl = dxy*(1/2+offset_yl); dw_yu = dxy*(1/2+offset_yu);
            for nn = 1:N_boundary
                ii = obj.inds_pipe_boundary(nn,1)+1;
                jj = obj.inds_pipe_boundary(nn,2)+1;
                N_faces = obj.VS_boundaries_2D(ii-1,jj-1); % Number of faces the cell has which comprise the interface
                
                iii_pipe = squeeze(G(ii,jj,2:end-1)); iii_left = squeeze(G(ii-1,jj,2:end-1)); iii_right = squeeze(G(ii+1,jj,2:end-1)); iii_bot = squeeze(G(ii,jj-1,2:end-1)); iii_top = squeeze(G(ii,jj+1,2:end-1));
                % At cell location
                K1 = k/dxy;
                q = q(end)+(1:Nz);
                iii(q) = iii_pipe;  jjj(q) = iii_pipe;  s(q) = K1; % At cell
                K3 = 0;
                
                % check left boundary
                if ~obj.VS_pipe_2D(ii-2,jj-1)
                    Bi = h./k*dw_xl;
                    % These values obtained from reformulating a boundary condition of the third kind into a*dT/dx+b*T=c
                    K2 = -k/dxy+h./(1+Bi);
                    K3 = K3 + h./(1+Bi).*(T_b-S.*(dw_xl)^2/2/k);
                else
                    K2 = 0;
                end
                q=q(end)+(1:Nz);
                iii(q) = iii_pipe;  jjj(q) = iii_left;  s(q) = K2/N_faces; % left of pipe.
                % check right boundary
                if ~obj.VS_pipe_2D(ii,jj-1)
                    Bi = h./k*dw_xu;
                    K2 = -k/dxy + h./(1+Bi); 
                    K3 = K3 + h./(1+Bi).*(T_b-S.*(dw_xu)^2/2/k);
                else
                    K2 = 0;
                end
                q=q(end)+(1:Nz);
                iii(q) = iii_pipe; jjj(q) = iii_right; s(q) = K2/N_faces; % right of pipe
                % check bottom boundary
                if ~obj.VS_pipe_2D(ii-1,jj-2)
                    Bi = h./k*dw_yl;
                    K2 = -k/dxy + h./(1+Bi);
                    K3 = K3 + h./(1+Bi).*(T_b-S.*(dw_yl)^2/2/k);
                else
                    K2 = 0;
                end
                q = q(end)+(1:Nz);
                iii(q) = iii_pipe;  jjj(q) = iii_bot;  s(q) = K2/N_faces; % bottom of pipe
            
                % check top boundary
                if ~obj.VS_pipe_2D(ii-1,jj)
                    Bi = h./k*dw_yu;
                    K2 = -k/dxy + h./(1+Bi);
                    K3 = K3 + h./(1+Bi).*(T_b-S.*(dw_yu)^2/2/k);
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
        %% getWallVariable
        function var_z = getWallVariable(obj,VS)
            TES_boundary_inds = find((6 - countConnectivity(~obj.VS_pipe_2D)).*(~obj.VS_pipe_2D));
            [ii,jj] = ind2sub(size(obj.VS_pipe_2D),TES_boundary_inds);
            N_TES_boundary = length(TES_boundary_inds);
            var_wall_all = zeros(obj.size_VS(3),N_TES_boundary);
            %ri = obj.radius;
            for nn = 1:N_TES_boundary
                var_wall_all(:,nn) = squeeze(VS(ii(nn),jj(nn),:)); % Column vector
                %ro = ((ii(nn)-0.5-obj.x_loc)^2 + ...
                 %   (jj(nn)-0.5-obj.y_loc)^2).^0.5;
               
            end
            var_z = mean(var_wall_all,2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% getWallTemperature

        function [T_wall_z] = getWallTemperature(obj,VS_T,thermal_conductivity,h_mean_z,T_fluid_z)
            xyw = obj.xy_width_vx; x_center = obj.x_loc; y_center = obj.y_loc; % I think this is a little more readable
            xl = x_center-xyw/2; xu = x_center+xyw/2; yl = y_center-xyw/2; yu = y_center+xyw/2;
            offset_xl = xl-round(xl); offset_xu = round(xu)-xu; offset_yl = yl-round(yl); offset_yu = round(yu)-yu;
            TES_boundary_inds = find((6 - countConnectivity(~obj.VS_pipe_2D)).*(~obj.VS_pipe_2D));
            [ii,jj] = ind2sub(size(obj.VS_pipe_2D),TES_boundary_inds);
            N_TES_boundary = length(TES_boundary_inds);
            boundary_x_ind_range = [floor(xl-0.5),ceil(xu+0.5)];
            boundary_y_ind_range = [floor(yl-0.5),ceil(yu+0.5)];
            T_wall_all = zeros(obj.size_VS(3),N_TES_boundary);
            dx = obj.vx_scale(1);
            dy = obj.vx_scale(2); % for now dx == dy but in case in the future i change that!\
            for nn = 1:N_TES_boundary
                T_cur = squeeze(VS_T(ii(nn),jj(nn),:)); % Column vector
                if ii(nn) == boundary_x_ind_range(1)
                    dw = (1/2+offset_xl)*dx;
                elseif ii(nn) == boundary_x_ind_range(2)
                    dw = (1/2+offset_xu)*dx;
                elseif jj(nn) == boundary_y_ind_range(1)
                    dw = (1/2+offset_yl)*dy;
                else
                    dw = (1/2+offset_yu)*dy;
                end
                if isscalar(thermal_conductivity)
                    k = thermal_conductivity;
                else
                    k = squeeze(thermal_conductivity(ii,jj,:));
                end
                Bi = dw*h_mean_z./k;
                %T_wall_all(:,nn) = T_cur - (T_cur-T_fluid_z)./(k./h_mean_z./(ri*dxy)+log(ro/ri))*log(ro/ri);
                T_wall_all(:,nn) = T_cur + Bi/(1+Bi)*(T_fluid_z-T_cur)+obj.S_wall_z.*dw^2./(2*k*(1+Bi));
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
                f_lam = 64/Re; 
            end
            if Re > 2100
                f_t = (-1.8*log10((rel_roughness/3.7).^1.11+6.9./Re)).^(-2); % Haaland correlation for turbulent region
            end
            if Re <= 2100
                f_Darcy = f_lam;
            elseif Re >= 4000
                f_Darcy = f_t;
            else
                alpha = (4000-Re)/(4000-2100); % Simple linear weighting of f_lam and f_t;
                f_Darcy = f_lam*alpha + (1-alpha)*f_t;
            end
        end
    end
    methods (Static)
        function v = getNuFitParams()
            v = [0.6995,0.6784,0.3606,0.5493,0.3116,53.3650];
        end
        function Nu_x = getNu_x_lam(x_star)
            % Correlation fit to data from Table 3.15, Bhatti and Shah 1987
            v = SquarePipe.getNuFitParams();
            a1 = v(1); b1 = v(2); c1 = v(3);
            a2 = v(4); b2 = v(5); c2 = v(6);
            if x_star <= 1/80 % 
                Nu_x = a1*x_star.^(-c1)+b1;
            else
                Nu_x = 2.975+a2.*x_star.^(-b2).*exp(-c2*x_star);
            end
        end
        function Nu_m = getNu_m_lam(x_star_L,x_star_U)
            % Correlation fit to data from Table 3.15, Bhatti and Shah
            % x_star = 1/Gz
            % x_star_L < x_star_U
            v = SquarePipe.getNuFitParams();
            a1 = v(1); b1 = v(2); c1 = v(3);
            a2 = v(4); b2 = v(5); c2 = v(6);
            if x_star_U <= 1/80
                Nu_m = integral(@(x) a1*x.^(-c1)+b1,x_star_L,x_star_U);
            elseif x_star_L >= 1/80
                Nu_m = integral(@(x) 2.975+a2.*x.^(-b2).*exp(-c2*x),1/1000,x_star_U);
            else % x_star_L < 1/80 && x_star_R > 1/80
                Nu_m = integral(@(x) a1*x.^(-c1)+b1,x_star_L,1/80) + integral(@(x) 2.975+a2*x.^(-b2).*exp(-c2*x),1/80,x_star_U);
            end
            Nu_m = Nu_m/(x_star_U-x_star_L);
        end
    end

end

