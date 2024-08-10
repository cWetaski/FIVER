clear
close all
clc

%% TES Parameters - borosilicate glass SHOTT-BK7

D_TES = 67/1000; % m;
L = 1.5*D_TES; % m
k_TES = 1.114; % W/(K-m); Shott BK7
c_TES = 858; % J/(kg K): specific heat
rho_TES = 2.51e3; % kg/m3
thermal_expansion = 8.3e-6; % /K
youngs_mod = 82e9; % Pa
poisson_ratio = 0.206;
tensile_strength = 832*9.81/100*1e6; % Pa conversion from kg/cm2 to Pa https://www.imetra.com/borosilicate-glass-material-properties/

%% Pipe parameters
roughness_pipe = 0.0015/1000; % m: for glass, source: https://www.engineersedge.com/fluid_flow/pipe-roughness.htm
T_air_i = 20 + 273.15; % Kelvin
m_dot = 0.002; % kg/s
N_grid_D = 16; % Number of gridcells across the diameter

%% Simulation parameters
T_TES_i = 400 + 273.15; % Kelvin: Initial temperature of the TES
t_sim = 60; % seconds
N_time_steps = 20;
time_step = t_sim/N_time_steps;

%% Derived parameters
alpha_TES = k_TES/(c_TES*rho_TES);

vx_scale = D_TES/N_grid_D; % side length of each grid cell
N_grid_L = round(L/vx_scale);

r_pipe = vx_scale*1;
r_vx_pipe = r_pipe/vx_scale;
D_pipe = r_pipe*2;
dz = vx_scale;
void_ratio = pi/4*D_pipe^2/(vx_scale^2); % percentage of cell that is empty due to pipe
V_pipe = void_ratio*L*vx_scale^2; % m3: total pipe volume

VS_T_init = zeros(N_grid_D,N_grid_D,N_grid_L);
size_VS = size(VS_T_init);

cylinder_bool = false(size_VS);
center = size_VS(2)/2;
for ii = 1:size_VS(1)
    x_dif = ii-0.5-center;
    for jj = 1:size_VS(2)
        y_dif = jj-0.5-center;
        R_sq = x_dif^2+y_dif^2;
        if R_sq <= (N_grid_D/2)^2
            cylinder_bool(ii,jj,:) = true;
        end
    end
end
VS_T_init(cylinder_bool) = T_TES_i; % Initial temperature field in TES, everything outside the TES can just be 0 for now

VS_V_TES = zeros(size_VS);
VS_V_TES(cylinder_bool) = vx_scale^3;

VS_T_fixed = false(size_VS);
VS_T_fixed(~cylinder_bool) = true;

VS_alpha = zeros(size_VS);
VS_alpha(cylinder_bool) = alpha_TES;

VS_opaq = false(size_VS);
VS_opaq_emissivities = VS_opaq;

voxel_space = VoxelSpace();

%% Define conduction problem;
mesh_3D = createMesh3D(size_VS(1),size_VS(2),size_VS(3),size_VS(1)*vx_scale,size_VS(2)*vx_scale,size_VS(3)*vx_scale);
BC = createBC(mesh_3D); % all Neumann (==0) boundary condition structure is default

alpha_cell = createCellVariable(mesh_3D, VS_alpha); % assign the thermal diffusivity to the cells
alpha_face = harmonicMean(alpha_cell); % calculate harmonic average of the diffusion coef on the cell faces
alpha_face.xvalue(isnan(alpha_face.xvalue))=0;% Get rid of nans that arise from adjacent cells having 0 diffusivity!
alpha_face.yvalue(isnan(alpha_face.yvalue))=0; 
alpha_face.zvalue(isnan(alpha_face.zvalue))=0;

% Get components of RHS and matrix of coefficients for diffusion and BCs (these do not change on iterations)
M_diff = diffusionTerm(alpha_face); % matrix of coefficients for the diffusion term
[M_BC, RHS_BC] = boundaryCondition(BC); % matrix of coefficients and RHS vector for the BC
M_size = size(M_diff);

%% Define our pipe(s) in the TES
N_pipe_side = 1;
N_pipe = N_pipe_side^2;
Nx = size_VS(1); Ny = size_VS(2); Nz = size_VS(3);
RHS_pipe = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
loc_pipes = zeros(N_pipe,2);
ind_pipes = zeros(N_pipe,2);
matrix_ind_map = reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2,Nz+2);

for n1 =1:N_pipe_side
    for n2 = 1:N_pipe_side
        nn = N_pipe_side*(n1-1)+n2;
        loc_pipes(nn,:) = round(2*[Nx*n1/(N_pipe_side+1),Ny*n2/(N_pipe_side+1)])/2;
        ind_pipes(nn,:) = floor(loc_pipes(nn,:))+1;
        if nn == 1
            plot_bool = true;
        else
            plot_bool = false;
        end
        [VS_pipe_2D{nn},VS_boundaries_2D{nn},inds_pipe_boundary{nn},inds_pipe_interior{nn},~] = getMultiPipe(size_VS,loc_pipes(nn,1),loc_pipes(nn,2),r_vx_pipe,plot_bool);
        VS_V_TES(repmat(VS_pipe_2D{nn},[1,1,size_VS(3)])) = 0; % don't include the pipe in the volume of TES
    end
end

RHS_pipe = zeros(M_size(1),1);
M_pipe = sparse(M_size(1),M_size(2));
T_w_z = T_TES_i;

for nn = 1:N_pipe
    inds = vertcat(inds_pipe_interior{nn},inds_pipe_boundary{nn})+1; % plus 1 to account for ghost cells
    for n2 = 1:size(inds,1)
        iii = squeeze(matrix_ind_map(inds(n2,1),inds(n2,2),2:end-1));
        for row_i = iii
            M_diff(iii,:) = 0; % Delete diffusion term from pipe locations
        end
    end
    [h_mean_z{nn},T_air_z{nn},Q_pipe_z,T_air_o] = GetQPipe(T_w_z,T_air_i,m_dot,vx_scale,D_pipe,roughness_pipe);
    [M_pipe_cur, RHS_pipe_cur] = PipeBC_3D_Multi(mesh_3D,loc_pipes(nn,1),loc_pipes(nn,2),r_vx_pipe,k_TES,h_mean_z{nn},T_air_z{nn});
    M_pipe = M_pipe + M_pipe_cur;
    RHS_pipe = RHS_pipe + RHS_pipe_cur;
end


%% Initialize counters
count_itr = 0; % Total iteration counter

%% define initial values
VS_T_old = VS_T_init;
VS_T_new = VS_T_old;
T_old = createCellVariable(mesh_3D, VS_T_init,BC); % initial values
source_cell = createCellVariable(mesh_3D,0); % Initialize source term cell variable

%% Set up fixed temperatures
% Pad VS_T_spec so it is the same size as the cell variables (which have ghost cells);
VS_T_spec_padded = padarray(VS_T_fixed,[1 1 1],0,'both');
fixed_inds = find(VS_T_spec_padded); % Get linear indices of fixed values
M_fixed = sparse(fixed_inds,fixed_inds,1, M_size(1), M_size(2)); % Create a sparse matrix of 0s except for diagonals where temperature is fixed
RHS_fixed = zeros(M_size(1),1);
RHS_fixed(fixed_inds) = T_old.value(fixed_inds); % Set fixed values on RHS

%% Results/Statistics
E_tot_i = sum(VS_T_old.*VS_V_TES*rho_TES*c_TES,'all'); % J: Initial energy content in TES

T_air_o = zeros(N_pipe,1);
T_air_o_t = zeros(N_time_steps,1);
Q_pipe_tot = zeros(N_pipe,1);
Q_pipe_tot_t = zeros(N_time_steps,1);
Delta_E_t = zeros(N_time_steps,1);
%M_diff_plus_pipe = M_diff;

%% Figure
T_2D_init_E = VS_T_init(:,:,1);
for nn = 1:N_pipe
    T_2D_init_E(VS_pipe_2D{nn}) = nan;
end
f = figure;
zrange = [(T_air_i-273.15),(T_TES_i-273.15)];
%xlim([0,Lx]); ylim([0,Ly]), zlim(zrange)
imagesc(T_2D_init_E-273.15);
colorbar();
clim(zrange);
%colormap jet
drawnow();

%% loop
for tt = 1:N_time_steps
    q_rel_diff = 1e6;
    count = 0;
    while q_rel_diff > 0.001 && count < 3
        count = count+1;
        [M_trans, RHS_trans] = transientTerm(T_old, time_step); % Get transient term based on old temperature field. Don't quite understand how this works but it works
        RHS_pipe = zeros(M_size(1),1);
        M_pipe = sparse(M_size(1),M_size(2));
        for n1 =1:N_pipe_side
            for n2 = 1:N_pipe_side
                nn = N_pipe_side*(n1-1)+n2;
                inds = vertcat(inds_pipe_interior{nn},inds_pipe_boundary{nn})+1; % plus 1 to account for ghost cells
                for n3 = 1:size(inds,1)
                    iii = squeeze(matrix_ind_map(inds(n3,1),inds(n3,2),2:end-1));
                    M_trans(sub2ind(size(M_trans),iii,iii)) = 0; % Delete transient term from pipe locations
                    RHS_trans(iii) = 0;
                end
                TES_boundary_inds = find((6 - countConnectivity(~VS_pipe_2D{nn})).*(~VS_pipe_2D{nn}));
                [ii_TES_boundary,jj_TES_boundary] = ind2sub(size_VS(1:2),TES_boundary_inds);
                N_TES_boundary = length(TES_boundary_inds);
                T_wall_all = zeros(size_VS(3),N_TES_boundary);
                for n3 = 1:N_TES_boundary
                    T_cur = squeeze(VS_T_old(ii_TES_boundary(n3),jj_TES_boundary(n3),:)); % Column vector
                    r_vx_cur = ((ii_TES_boundary(n3)-0.5-loc_pipes(nn,1))^2 + ...
                             (jj_TES_boundary(n3)-0.5-loc_pipes(nn,2))^2).^0.5;
                    T_wall_all(:,n3) = T_cur - (T_cur-T_air_z{nn}(:))./(k_TES./h_mean_z{nn}(:)./r_pipe+log(r_vx_cur/r_vx_pipe))*log(r_vx_cur/r_vx_pipe);
                    T_boundary_all(:,n3) = T_cur;
                    r_boundary_all(1,n3) = r_vx_cur; 
                end
                r_boundary = mean(r_boundary_all);
                weights = 1 - (r_boundary_all-r_boundary)^2/r_boundary^2;
                T_boundary = sum(T_boundary_all.*weights,2)/sum(weights);
                T_w_z = mean(T_wall_all,2);

                [thermal_stresses,vm_stresses] = getThermalStressCylinder(r_vx_pipe,r_vx_pipe,r_boundary,T_w_z-T_boundary,thermal_expansion,poisson_ratio,youngs_mod);
                max(vm_stresses/tensile_strength)
                [h_mean_z{nn},T_air_z{nn},Q_pipe_z,T_air_o(nn)] = GetQPipe(T_w_z,T_air_i,m_dot,vx_scale,D_pipe,roughness_pipe);

                Q_pipe_tot(nn) = sum(Q_pipe_z);
                [M_pipe_cur, RHS_pipe_cur] = PipeBC_3D_Multi(mesh_3D,loc_pipes(nn,1),loc_pipes(nn,2),r_vx_pipe,k_TES,h_mean_z{nn},T_air_z{nn});    
                M_pipe = M_pipe+M_pipe_cur;
                RHS_pipe = RHS_pipe + RHS_pipe_cur;
            end
        end
        air_table = readtable('air_properties.txt','Delimiter',' '); % Table of air properties vs T at p=1atm

        T_air_o_t(tt) = mean(T_air_o);
        [~,cp_mean,~] = GetAirProps((T_air_o_t(tt)+T_air_i)/2,air_table); % Get initial cp

        Q_air_t(tt) = -(T_air_o_t(tt)-T_air_i)*cp_mean*m_dot*N_pipe;
        Q_pipe_tot_t(tt) = -sum(Q_pipe_tot);

        M = M_trans-M_diff+M_pipe+M_BC;
        RHS = RHS_trans+RHS_BC+RHS_pipe;
        
        % Set fixed temperature values;
        M(fixed_inds,:) = 0;
        M = M + M_fixed;
    
        RHS(fixed_inds) = 0;
        RHS = RHS + RHS_fixed;
    
        T_new = solvePDE(mesh_3D,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS but returns a cell variable
        VS_T_new = T_new.value(2:(end-1),2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)
        Delta_E_t(tt) = sum((VS_T_new-VS_T_old).*VS_V_TES*rho_TES*c_TES,'all'); % Joules
        Q_TES_t(tt) = Delta_E_t(tt)/time_step;
        q_rel_diff = abs(Q_TES_t(tt)-Q_pipe_tot_t(tt))/abs(Q_pipe_tot_t(tt));
    end
    %count
    VS_T_old = VS_T_new; % Update temperature field for internal iteration (or time step, if internal iteration is finished
    % Update temperature field for FVTool transient term
    T_old = T_new;
    % Update plot
    T_2D_new_E = VS_T_new(:,:,1); % Get a 2D cross section
    for nn = 1:N_pipe
        T_2D_new_E(VS_pipe_2D{nn}) = nan;
    end
    if true%mod(tt/N_time_steps*100,1) == 0
        %T_2D_wall_new = add_wall_T_2D(T_2D_new,loc_pipes);
        %surf(xx_vals,yy_vals,T_2D_wall_new-273.15,edgeColor='none');
        %view(0,90);
        %xlim([0,Lx]); ylim([0,Ly]); zlim(zrange);
        %text(Lx*0.75,zrange(1)+(zrange(2)-zrange(1))*0.75,sprintf('t = %0.1f',tt*time_step));
        imagesc(T_2D_new_E-273.15);
        colorbar()
        clim(zrange);
        title(sprintf("Time = %0.1f s",time_step*tt))
        drawnow();
        pause(0.01)
        %T_2D_wall_old = T_2D_wall_new;
        %fprintf("time = %0.1f\n",tt*time_step)
    end
    
end

E_tot_f = sum(VS_T_new.*VS_V_TES*rho_TES*c_TES,'all'); % J: Final energy content in TES

fprintf('Energy change in TES: %0.4f kJ \n',(E_tot_f-E_tot_i)/1000);
fprintf('Average heat transfer rate: %0.4f W \n',mean(Q_pipe_tot_t));
fprintf('Energy change according to average heat transfer rate: %0.4f kJ \n',t_sim*mean(Q_pipe_tot_t)/1000);
t_vals = (1:N_time_steps)*time_step;

f = figure;
hold on
xlabel("Time (s)")
yyaxis left
plot(t_vals,T_air_o_t-273.15);
ylabel('Outlet temperature')
yyaxis right
plot(t_vals,Q_pipe_tot_t,'--');
plot(t_vals,Delta_E_t/time_step,'-')
plot(t_vals,Q_air_t,':')
ylabel('Heat transfer rate (W)')

fprintf('%0.4f\n',Delta_E_t(1)/time_step-(Q_pipe_tot_t(1)))

f = figure;
plot(t_vals,Q_pipe_tot_t-Delta_E_t/time_step)


% 1000 seconds over 1000 time steps with internal iteration (0.001) == -113.1356 kJ (TES)  [Note: internal iteration only for first few iterations] 
% 1000 seconds over 50 time steps with internal iteration (0.001) == -113.0409 kJ (TES)
% 1000 seconds over 50 time steps without internal iteration == -112.7774 kJ

% 2 hours over 50 time steps without internal iteration == -677.6357 kJ
% 2 hours over 50 time steps with internal iteration (0.001) == -679.5247 kJ
% 2 hours over 7000 time steps == -680.1024 kJ

% 100 seconds over 1000 time steps without internal iteration == -12.8841 kJ
% 100 seconds over 1000 time steps with internal iteration (0.001) == -12.8850 kJ
% 100 seconds over 5 time steps without internal iteration == -12.6013 kJ
% 100 seconds over 5 time steps with internal iteration (0.001) == -12.8240 kJ
