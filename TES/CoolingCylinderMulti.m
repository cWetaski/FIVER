clear
close all
clc

Bi_cyl = 1;
Fo_slab_step = 1;
N_stop = 1;
Nxy_cylinder = 42; % Voxel diameter of the cylinder

ro = 2; % m: cylinder inner diameter
Do = ro*2; % m: (cylinder outer diameter)
dxy = Do/Nxy_cylinder; % voxel length
ro_vx = ro/dxy;

ri_vx = 1; % inner diameter
ri = ri_vx*dxy;

q_gen = 0.04; % W/m3: Volumetric heating rate scale in this way with voxel diameter to have roughly invariant temperature when changing resolution
h_mean = 1; % W/(K-m2)
k_cyl = h_mean*ri/Bi_cyl; % W/(K-m)
k_cyl = 1;

alpha_t = Fo_slab_step*ro^2; %m/s
time_step = 100; % s
alpha_cyl = alpha_t/time_step; % m2/s
alpha_cyl = 1;
N_time_steps = 200;
t_sim = time_step*N_time_steps;

Bi_cell = Bi_cyl*dxy/ri;
Fo_cell_step = Fo_slab_step*(ro/dxy)^2;

rho_cyl = 1; % kg/m3
c_cyl = k_cyl/rho_cyl/alpha_cyl; % J/kg-K
T_air = 0; % K
Bi_cyl = h_mean*ri/Bi_cyl;

%% Analytic solution
r_vx = 1:round(Nxy_cylinder/2);
r_vx((r_vx-0.5)<ri_vx) = [];
r_vals = r_vx*dxy;
r_vals(r_vals>ro) = [];
R_vals = r_vals/ri;
N_rr = length(R_vals);

% Analytic temp
N_analytic = 100;
Ro = ro/ri;
R_vals_analytic = linspace(1,Ro,N_analytic);
sol = HollowCylinderHeatGen(R_vals_analytic,Ro,Bi_cyl);
T_R_analytic = T_air + sol*q_gen*ri^2/k_cyl;
T_max = T_R_analytic(end);
T_max = 1;
T_R_analytic = T_R_analytic/T_max; % Scale to interval [0,1]
T_cyl_i = T_max;


%% 2D conduction problem with interior convective BC
Nx = Nxy_cylinder; Ny = Nxy_cylinder;
size_VS = [Nx,Ny];
Lx = Do; Ly = Do;
dxy = Lx/Nx;

cylinder_bool = false(Nx,Ny);
center = [Nx,Ny]/2;
max_R = 0;
for ii = 1:Nx
    x_dif = ii-0.5-center(1);
    for jj = 1:Ny
        y_dif = jj-0.5-center(2);
        rr_sq = x_dif^2+y_dif^2;
        
        if rr_sq <= ro_vx^2
            if sqrt(rr_sq) > max_R
                max_R = sqrt(rr_sq);
            end
            cylinder_bool(ii,jj) = true;
        end
    end
end

mesh_2D = createMesh2D(Nx,Ny,Lx,Ly);
BC_2D = createBC(mesh_2D);
BC_2D.left.a = 1; BC_2D.left.b = 0; BC_2D.left.c=0;
BC_2D.right.a = 1; BC_2D.right.b = 0; BC_2D.right.c = 0;
BC_2D.bottom.a = 1; BC_2D.bottom.b = 0; BC_2D.bottom.c = 0;
BC_2D.top.a = 1; BC_2D.top.b = 0; BC_2D.top.c = 0;
[M_BC,RHS_BC] = boundaryCondition(BC_2D);

%% Diffusion term
alpha_2D = zeros(Nx,Ny);
alpha_2D(cylinder_bool) = alpha_cyl;
alpha_cell_2D = createCellVariable(mesh_2D,alpha_2D);
alpha_face_2D = harmonicMean(alpha_cell_2D);
alpha_face_2D.xvalue(isnan(alpha_face_2D.xvalue))=0;% Get rid of nans that arise from adjacent cells having 0 diffusivity!
alpha_face_2D.yvalue(isnan(alpha_face_2D.yvalue))=0; 
M_diff = diffusionTerm(alpha_face_2D);

%% Apply pipe to system
M_pipe = M_diff;
RHS_pipe = zeros((Nx+2)*(Ny+2),1);

loc_pipe = [Nx/2,Ny/2];
ind_pipe = floor(loc_pipe)+1;
matrix_ind_map = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
    
[VS_pipe,VS_boundaries,inds_pipe_boundary,inds_pipe_interior,~] = getMultiPipe(size_VS,loc_pipe(1),loc_pipe(2),ri_vx,true);
[M_pipe, RHS_pipe] = PipeBC_2D_Multi(mesh_2D,loc_pipe(1),loc_pipe(1),ri_vx,k_cyl,h_mean,T_air);

inds = vertcat(inds_pipe_boundary,inds_pipe_interior);
for nn = 1:size(inds,1)
    M_diff(matrix_ind_map(inds(nn,1)+1,inds(nn,2)+1),:) = 0; % Delete diffusion term from pipe locations
end


%% Apply source term to system
source_term_2D = zeros(Nx,Ny);
source_term_2D(cylinder_bool) = q_gen/rho_cyl/c_cyl; % T/s
source_term_2D(VS_pipe)=0;

source_cell_2D = createCellVariable(mesh_2D,source_term_2D);
RHS_source = constantSourceTerm(source_cell_2D);

%% initial Temperature values
T_2D_init = zeros(Nx,Ny);
T_2D_init(cylinder_bool) = T_cyl_i;
T_2D_cell_init = createCellVariable(mesh_2D, T_2D_init,BC_2D); % initial values

T_2D_cell_old = T_2D_cell_init;
T_2D_init_E = T_2D_init;
T_2D_init_E(sub2ind([Nx,Ny],ind_pipe(:,1),ind_pipe(:,2))) = 0; % We don't want the ghost cell value to contribute to energy calculation!
T_2D_old_E = T_2D_init_E;

x_vals = (1:Nx)*dxy-dxy*0.5;
y_vals = (1:Nx)*dxy-dxy*0.5;
count = 0;

inds = vertcat(inds_pipe_interior,inds_pipe_boundary);
VS_T_new = T_2D_cell_old.value(2:end-1,2:end-1);
VS_T_new(sub2ind(size(VS_T_new),inds(:,1),inds(:,2))) = nan;

f = figure();
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact'); 
f.Position =  [133.8000  242.0000  914.2000  420.0000];
nexttile();
plot(R_vals_analytic,T_R_analytic,'--','Color','k');
T_R = ones(N_rr,1);
xrange = [0,max(R_vals)*1.1];
zrange = [T_air,1.1];
hold on;
plot(R_vals,T_R);
xline(1,':','R_i = 1','Color','k')
xline(Ro,':',sprintf('R_o = %0.1f',Ro),'Color','k');
ylim(zrange);
xlim(xrange);
%xlim([0,Lx]); ylim([0,Ly]), zlim(zrange)
%imagesc(T_2D_init_E);

%colormap jet
drawnow();
hold off;
Q_wall = h_mean*(T_air-T_cyl_i)*pi*2*r_vx*dxy;
matrix_ind_map = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

for tt = 1:N_time_steps
    VS_T_old = VS_T_new;
    [M_trans, RHS_trans] = transientTerm(T_2D_cell_old, time_step); % Get transient term based on old temperature field. Don't quite understand how this works but it works

    inds = vertcat(inds_pipe_interior,inds_pipe_boundary);
    ii = inds(:,1)+1;
    jj = inds(:,2)+1;
    for n2 = 1:length(ii)
        iii(n2,1) = matrix_ind_map(ii(n2),jj(n2));
    end
    jjj = iii;
    RHS_trans(iii) = 0;
    M_trans(iii,jjj) = 0;
    T = T_2D_cell_old.value;
    %T_wall(n1) = 1/8*(T(ii,jj)*4 + T(ii-1,jj) + T(ii+1,jj)+T(ii,jj-1)+T(ii,jj+1));
    %Q_wall(n1) = h_mean*(T_air-T_wall(n1))*pi*2*r_vx*dxy; % W/m: Negative value indicates heat leaving the TES, multiply by 4 since we have 4 walls for each pipe, multiply by pi/4 since the pipe is cylindrical, not flat wall

    %q_wall_tot(tt) = sum(q_wall);
 
    M = M_trans-M_diff+M_pipe + M_BC;
    RHS = RHS_trans+RHS_BC + RHS_pipe + RHS_source;

    T_2D_cell_new = solvePDE(mesh_2D,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS
    VS_T_new = T_2D_cell_new.value(2:end-1,2:end-1);
    T_ghost = VS_T_new(sub2ind(size(VS_T_new),inds(:,1),inds(:,2)));
    T_ghost_mean = mean(T_ghost);
    VS_T_new(sub2ind(size(VS_T_new),inds(:,1),inds(:,2))) = nan;

    % Update temperature field for FVTool transient term
    nexttile(2);
    imagesc(VS_T_new);
    colorbar();
    T_2D_cell_old = T_2D_cell_new;
    Q = sum(VS_T_old(:)-VS_T_new(:),'omitnan')*c_cyl*rho_cyl*(sum(cylinder_bool(:))-sum(VS_pipe(:)))*dxy^2/time_step;
    VS_cyl = ~VS_pipe;
    VS_cyl_boundary = 6-countConnectivity(VS_cyl);
    VS_cyl_boundary(VS_cyl_boundary==6) = 0;
    VS_cyl_boundary = VS_cyl_boundary>0; % convert to logical
    T_boundary_old = mean(VS_T_old(VS_cyl_boundary),'all');
    ind_cyl_boundary = find(VS_cyl_boundary);
    [ii_cyl_boundary,jj_cyl_boundary] = ind2sub(size(VS_T_new),ind_cyl_boundary);
    for nn = 1:length(ii_cyl_boundary)
        r_cyl_boundary(nn) = ((ii_cyl_boundary(nn)-0.5-center(1))^2+(jj_cyl_boundary(nn)-0.5-center(2))^2)^0.5;
    end
    ro_boundary = mean(r_cyl_boundary)*dxy;
    Tw_old = T_boundary_old - (T_boundary_old-T_air)/(k_cyl/h_mean/ri+log(ro_boundary/ri))*log((ro_boundary/ri));
    Q2 = h_mean*(Tw_old-T_air)*(2*ri*pi)^2/(sum(VS_boundaries(:))*dxy)
    Q3 = k_cyl*(T_boundary_old - T_ghost_mean)/dxy
 
    
    T_2D_new = T_2D_cell_new.value(2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)
    T_R = zeros(N_rr,1);
    counts = zeros(N_rr,1);
    for ii = 1:Nx
        x_dif = ii-0.5-center(1);
        for jj = 1:Ny
            y_dif = jj-0.5-center(2);
            if VS_pipe(ii,jj)
                continue
            end
            if ii == 20 && jj == 21
                debug =0;
            end
            d_r = sqrt(x_dif^2+y_dif^2);
            if d_r > ro_vx || d_r < ri_vx
                continue
            end
            ind = round(d_r-ri_vx);
            T_R(ind) = T_R(ind) + T_2D_new(ii,jj);
            counts(ind) = counts(ind)+1;
        end
    end
    T_R = T_R./counts/T_max;        

    % Change in energy
    T_2D_new_E = T_2D_new;
    T_2D_new_E(VS_pipe) = 0; % We don't want the ghost cell value to contribute to energy calculation!
    Delta_E(tt) = sum(T_2D_new_E-T_2D_old_E,'all').*dxy^3*rho_cyl*c_cyl; % Joules/m
    T_2D_old_E = T_2D_new_E;

    % Update plot
    if mod(tt,N_stop)==0%mod(tt/N_time_steps*100,1) == 0
        %T_2D_wall_new = add_wall_T_2D(T_2D_new,loc_pipes);
        %surf(xx_vals,yy_vals,T_2D_wall_new-273.15,edgeColor='none');
        %view(0,90);
        %xlim([0,Lx]); ylim([0,Ly]); zlim(zrange);
        %text(Lx*0.75,zrange(1)+(zrange(2)-zrange(1))*0.75,sprintf('t = %0.1f',tt*time_step));
        %imagesc(T_2D_new_E);
        %colorbar()
        %clim(zrange);
        nexttile(1);
        plot(R_vals_analytic,T_R_analytic,'--','Color','k');

        ylim(zrange);
        xlim(xrange);
        hold on
        plot(R_vals,T_R);
        xline(1,':','R_i = 1','Color','k')
        xline(Ro,':',sprintf('R_o = %d',Ro),'Color','k');
        legend({'Analytic','CFD'},'Location','southeast')
        xlabel('R = r/r_i')
        ylabel('T/ T_{o,analytic}')
        drawnow()
        hold off
        pause(0.01)
        %T_2D_wall_old = T_2D_wall_new;
    end
end

Delta_E_tot = (sum(T_2D_new_E-T_2D_init_E,'all')*dxy^3*rho_cyl*c_cyl); % Joules/m
%Delta_E_from_q = mean(q_wall_tot)*t_sim
%Delta_E_from_q = mean(q_wall_tot(1:end-1)+q_wall_tot(2:end))/2*t_sim

function sol = HollowCylinderHeatGen(R,Ro,Bi)
    sol = 1/4*(2/Bi.*(Ro.^2-1)+1-R.^2+2*Ro.^2.*log(R));
end