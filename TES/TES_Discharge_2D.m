clear
close all
clc

%% TES Parameters
k_TES = 5.4; % W/(K-m);https://technicalglass.com/technical-properties/
c_TES = 670; % J/(kg K): % specific heat
rho_TES = 2.2e3; % kg/m3;

%% Simulation parameters
T_TES_i = 200 + 273.15; % Kelvin: Initial temperature of the TES
t_sim = 100; % seconds
N_time_steps = 20;
time_step = t_sim/N_time_steps;

%% Derived parameters
alpha_TES = k_TES/(c_TES*rho_TES);

%% 2D conduction problem with interior convective BC
Nx = 25;
Ny = Nx;
Lx = 0.5;
Ly = Lx;
dxy = Lx/Nx;
h_mean = 100; % W/m2-K
T_air = 50+273.15; % K

mesh_2D = createMesh2D(Nx,Ny,Lx,Ly);
BC_2D = createBC(mesh_2D);
BC_2D.left.a = 1; BC_2D.left.b = 0; BC_2D.left.c=0;
BC_2D.right.a = 1; BC_2D.right.b = 0; BC_2D.right.c = 0;
BC_2D.bottom.a = 1; BC_2D.bottom.b = 0; BC_2D.bottom.c = 0;
BC_2D.top.a = 1; BC_2D.top.b = 0; BC_2D.top.c = 0;
[M_BC,RHS_BC] = boundaryCondition(BC_2D);

alpha_2D = alpha_TES;
alpha_cell_2D = createCellVariable(mesh_2D,alpha_2D);
alpha_face_2D = harmonicMean(alpha_cell_2D);
M_diff = diffusionTerm(alpha_face_2D);

%% initial values
T_2D = ones(Nx,Ny)*T_TES_i;
T_2D_cell_init = createCellVariable(mesh_2D, T_2D,BC_2D); % initial values

%% Apply pipe to system
N_pipe_side = 1;
N_pipe = N_pipe_side^2;
T_2D_init = T_2D;
size_M = size(M_diff);
RHS_pipe = zeros(size_M(1),1);
M_pipe = sparse(size_M(1),size_M(2));
loc_pipes = zeros(N_pipe,2);

matrix_ind_map = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
for n1 =1:N_pipe_side
    for n2 = 1:N_pipe_side
        nn = N_pipe_side*(n1-1)+n2;
        loc_pipes(nn,:) = [round(Nx*n1/(N_pipe_side+1)+0.5),round(Ny*n2/(N_pipe_side+1)+0.5)];
        ii = loc_pipes(nn,1);
        jj = loc_pipes(nn,2);
        [M_pipe_cur, RHS_pipe_cur] = PipeBC_2D(mesh_2D,ii,jj,k_TES,h_mean,T_air);
        M_pipe = M_pipe + M_pipe_cur;
        RHS_pipe = RHS_pipe + RHS_pipe_cur;

        iii = matrix_ind_map(ii+1,jj+1);
        M_diff(iii,:) = 0;
    end
end
T_2D_cell_old = T_2D_cell_init;
T_2D_init_E = T_2D_init;
T_2D_init_E(sub2ind([Nx,Ny],loc_pipes(:,1),loc_pipes(:,2))) = 0; % We don't want the ghost cell value to contribute to energy calculation!
T_2D_old_E = T_2D_init_E;

f = figure;
zrange = [(T_air-273.15),(T_TES_i-273.15)];
imagesc(T_2D_init_E-273.15);
drawnow();
q_wall = 4*h_mean*pi/4*(T_air-T_TES_i);
for tt = 1:N_time_steps
    [M_trans, RHS_trans] = transientTerm(T_2D_cell_old, time_step); % Get transient term based on old temperature field. Don't quite understand how this works but it works
    for n1 = 1:N_pipe
        ii = loc_pipes(n1,1)+1;
        jj = loc_pipes(n1,1)+1;
        iii = sub2ind([Nx+2,Ny+2],ii,jj);
        jjj = iii;
        RHS_trans(iii) = 0;
        M_trans(iii,jjj) = 0;
        T = T_2D_cell_old.value;
        T_wall(n1) = 1/8*(T(ii,jj)*4 + T(ii-1,jj) + T(ii+1,jj)+T(ii,jj-1)+T(ii,jj+1));
        q_wall(n1) = 4*h_mean*pi/4*(T_air-T_wall(n1))*dxy^2; % W/m: Negative value indicates heat leaving the TES, multiply by 4 since we have 4 walls for each pipe, multiply by pi/4 since the pipe is cylindrical, not flat wall
    end
    q_wall_tot(tt) = sum(q_wall);
    
  

    M = M_trans-M_diff+M_pipe + M_BC;
    RHS = RHS_trans+RHS_BC + RHS_pipe;

    T_2D_cell_new = solvePDE(mesh_2D,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS
    % Update temperature field for FVTool transient term
    T_2D_cell_old = T_2D_cell_new;
    
    T_2D_new = T_2D_cell_new.value(2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)
    
    % Change in energy
    T_2D_new_E = T_2D_new;
    T_2D_new_E(sub2ind([Nx,Ny],loc_pipes(:,1),loc_pipes(:,2))) = 0; % We don't want the ghost cell value to contribute to energy calculation!
    Delta_E(tt) = sum(T_2D_new_E-T_2D_old_E,'all').*dxy^3*rho_TES*c_TES; % Joules/m
    T_2D_old_E = T_2D_new_E;

    % Update plot
    if true%mod(tt/N_time_steps*100,1) == 0
        %T_2D_wall_new = add_wall_T_2D(T_2D_new,loc_pipes);
        %surf(xx_vals,yy_vals,T_2D_wall_new-273.15,edgeColor='none');
        %view(0,90);
        %xlim([0,Lx]); ylim([0,Ly]); zlim(zrange);
        %text(Lx*0.75,zrange(1)+(zrange(2)-zrange(1))*0.75,sprintf('t = %0.1f',tt*time_step));
        imagesc(T_2D_new_E-273.15);
        colorbar()
        clim(zrange);
        drawnow();
        pause(0.01)
        %T_2D_wall_old = T_2D_wall_new;
    end
end

Delta_E_tot = (sum(T_2D_new_E-T_2D_init_E,'all')*dxy^3*rho_TES*c_TES) % Joules/m
Delta_E_from_q = mean(q_wall_tot)*t_sim

Delta_E_from_q = mean(q_wall_tot(1:end-1)+q_wall_tot(2:end))/2*t_sim