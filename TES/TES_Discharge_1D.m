clear
close all
clc

%% TES Parameters
k_TES = 1.4; % W/(K-m);https://technicalglass.com/technical-properties/
c_TES = 670; % J/(kg K): % specific heat
rho_TES = 2.2e3; % kg/m3;

%% Simulation parameters
T_TES_i = 200 + 273.15; % Kelvin: Initial temperature of the TES
t_sim = 100; % seconds
N_time_steps = 100;
time_step = t_sim/N_time_steps;

%% Derived parameters
alpha_TES = k_TES/(c_TES*rho_TES);

%% 1D conduction problem with interior convective BC
Nx = 11;
Lx = 0.01;
dx = Lx/Nx;
h_mean = 45; % W/m2-K
T_air = 50+273.15; % K

mesh_1D = createMesh1D(Nx,Lx);
BC_1D = createBC(mesh_1D);
BC_1D.left.a = 1; BC_1D.left.b = 0; BC_1D.left.c=0;
BC_1D.right.a = 1; BC_1D.right.b = 0; BC_1D.right.c = 0;
[M_BC,RHS_BC] = boundaryCondition(BC_1D);

alpha_1D = max(alpha_TES);
alpha_cell_1D = createCellVariable(mesh_1D,alpha_1D);
alpha_face_1D = harmonicMean(alpha_cell_1D);
M_diff = diffusionTerm(alpha_face_1D);

%% initial values
T_1D = ones(Nx,1)*T_TES_i;
T_1D_cell_init = createCellVariable(mesh_1D, T_1D,BC_1D); % initial values


%% Apply pipe to system
N_pipe = 1;
T_1D_init = T_1D;
size_M = size(M_diff);
RHS_pipe = zeros(size_M(1),1);
M_pipe = sparse(size_M(1),size_M(2));
matrix_ind_map = 1:size_M(1);
for n =1:N_pipe
    ii_pipes(n) = round(Nx*n/(N_pipe+1)+0.5);
    [M_pipe_cur, RHS_pipe_cur] = PipeBC_1D(mesh_1D,ii_pipes(n),k_TES,h_mean,T_air);
    M_pipe = M_pipe + M_pipe_cur;
    RHS_pipe = RHS_pipe + RHS_pipe_cur;
    iii_pipe = matrix_ind_map(ii_pipes(n)+1);
    M_diff(iii_pipe,:)= 0;
end
T_1D_cell_old = T_1D_cell_init;
T_1D_old = T_1D_cell_init.value;
T_1D_init_E = T_1D_init;
T_1D_init_E(ii_pipes) = 0; % We don't want the ghost cell value to contribute to energy calculation!
T_1D_old_E = T_1D_init_E;

x_vals = (1:Nx)*dx-dx*0.5;
count = 0;
for ii_pipe = ii_pipes
    ii = ii_pipe+count;
    x_vals = [x_vals(1:(ii-1)),x_vals(ii)-dx/2,x_vals(ii)+dx/2,x_vals((ii+1):end)];
    count = count+1;
end

f = figure;
T_1D_wall_init = add_wall_T(T_1D_init(:),ii_pipes)'; % For plotting only!
T_1D_wall_old = T_1D_wall_init;
plot(x_vals,T_1D_wall_init-273.15);
yrange = [(T_air-273.15)*0.95,(T_TES_i-273.15)*1.05];
xlim([Lx*0.4,Lx*0.6]); ylim(yrange)
drawnow();
q_wall = 2*h_mean*pi/4*(T_air-T_TES_i);
for tt = 1:N_time_steps
    [M_trans, RHS_trans] = transientTerm(T_1D_cell_old, time_step); % Get transient term based on old temperature field. Don't quite understand how this works but it works
    count = 0;
    for ii_pipe = ii_pipes
        RHS_trans(ii_pipe+1) = 0;
        M_trans(ii_pipe+1,ii_pipe+1) = 0;

        count = count+1;
    end
    T_wall = 1/(4*N_pipe)*(T_1D_cell_old.value(ii_pipes) + T_1D_cell_old.value(ii_pipes+1)*2 + T_1D_cell_old.value(ii_pipes));
    q_wall = 2*h_mean*pi/4*(T_air-T_wall)*dx; % w/m: Negative value indicates heat leaving the TES, multiply by 2 since we have two walls for each pipe, multiply by pi/4 since the pipe is cylindrical, not flat wall
    q_wall_tot(tt+1) = sum(q_wall);
  

    M = M_trans-M_diff+M_pipe+M_BC; %+ M_pipe_BC;
    RHS = RHS_trans+RHS_BC + RHS_pipe;

    T_1D_cell_new = solvePDE(mesh_1D,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS
    % Update temperature field for FVTool transient term
    T_1D_cell_old = T_1D_cell_new;
    
    T_1D_new = T_1D_cell_new.value(2:(end-1)); % updating temperatures (ignoring ghost cells)
    
    

    
    % Change in energy
    T_1D_new_E = T_1D_new;
    T_1D_new_E(ii_pipes) = 0;
    Delta_E(tt) = sum((T_1D_new_E-T_1D_old_E).*dx^2*rho_TES*c_TES); % Joules/m
    T_1D_old_E = T_1D_new_E;

    % Update plot
    if mod(tt/N_time_steps*100,1) == 0
        T_1D_wall_new = add_wall_T(T_1D_new(:),ii_pipes)';
        plot(x_vals,T_1D_wall_new-273.15);
        xlim([0,Lx]); ylim(yrange)
        text(Lx*0.75,yrange(1)+(yrange(2)-yrange(1))*0.75,sprintf('t = %0.1f',tt*time_step));
        drawnow();
        pause(0.05)
        T_1D_wall_old = T_1D_wall_new;
    end
end

Delta_E_tot = sum((T_1D_new_E-T_1D_init_E).*dx^2*rho_TES*c_TES) % Joules/m
Delta_E_from_q = mean(q_wall_tot(1:end-1)+q_wall_tot(2:end))/2*t_sim

function T = add_wall_T(T,ii_pipes)
    count = 0;
    for ii_pipe = ii_pipes
        ii = ii_pipe+count;
        T = [T(1:(ii-1))',(T(ii-1)+T(ii))/2,(T(ii+1)+T(ii))/2,T((ii+1):end)']';
        count = count+1;
    end
end
