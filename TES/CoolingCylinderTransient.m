clear
close all
clc

Bi_i = 1;
Bi_o = 0;
Fo_step = 0.01;
N_stop = 1;
Nxy_cylinder = 51; % Voxel diameter of the cylinder
save_animation = true;

ro = 2; % m: cylinder inner diameter
Do = ro*2; % m: (cylinder outer diameter)
dxy = Do/Nxy_cylinder; % voxel length

ri_vx = 0.5;
ri = dxy*ri_vx; % m: Inner diameter is (currently) only 1 voxel wide

ro_vx = ro/dxy;

h_mean = 100; % W/(K-m2)
k_cyl = h_mean*ri/Bi_i; % W/(K-m)

alpha_t = Fo_step*ro^2; %m/s
time_step = 100; % s
alpha_cyl = alpha_t/time_step; % m2/s
N_time_steps = 400;
t_sim = time_step*N_time_steps;

Bi_cell = Bi_i*dxy/ri;
Fo_cell_step = Fo_step*(ro/dxy)^2;

rho_cyl = 1000; % kg/m3
c_cyl = k_cyl/rho_cyl/alpha_cyl; % J/kg-K
T_air = 0; % K
T_cyl_0 = 1; % K: Initial temperature

%% Analytic solution
N_rr = round(Nxy_cylinder-1)/2;
R_vals = ((1:N_rr))*2;
r_vals = R_vals*ri;
r_vals_analytic = [dxy/2,r_vals,r_vals(end)+dxy/2];

% Analytic temp
R_vals_analytic = r_vals_analytic/ri;
Ro = ro/ri;
N_lambda = 40;

lambdas = get_lambdas(N_lambda,Bi_i,Bi_o,Ro);
W = get_W(R_vals_analytic,lambdas,Bi_i);
E = get_E(lambdas,Bi_i,Bi_o,T_cyl_0,T_air,Ro);
sol = HollowCylinderConvectionTransient(R_vals_analytic,0,W,E);
T_R_analytic = T_air + sol;

%% 2D conduction problem with interior convective BC
Nx = Nxy_cylinder; Ny = Nxy_cylinder;
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
        
        if rr_sq <= (Nx/2)^2
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
loc_pipe = [Nx/2,Ny/2];
ind_pipe = floor(loc_pipe)+1;
matrix_ind_map = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

M_diff(matrix_ind_map(ind_pipe(1)+1,ind_pipe(2)+1),:) = 0; % Removing diffusion from row associated with the pipe           
[M_pipe, RHS_pipe] = PipeBC_2D(mesh_2D,ind_pipe(1),ind_pipe(2),k_cyl,h_mean,T_air);


%% initial Temperature values
T_2D_init = zeros(Nx,Ny);
T_2D_init(cylinder_bool) = T_cyl_0;
T_2D_cell_init = createCellVariable(mesh_2D, T_2D_init,BC_2D); % initial values

T_2D_cell_old = T_2D_cell_init;
T_2D_init_E = T_2D_init;
T_2D_init_E(sub2ind([Nx,Ny],ind_pipe(1),ind_pipe(2))) = 0; % We don't want the ghost cell value to contribute to energy calculation!
T_2D_old_E = T_2D_init_E;

if save_animation
    myVideo = VideoWriter('CoolingCylinderTransientAnim','MPEG-4'); %open video file
    myVideo.FrameRate = 24;  %can adjust this, 5 - 10 works well for me
    open(myVideo);
end

f = figure;
colormap(f,getPyPlot_cMap("inferno"));
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
f.Position = [133.8000  242.0000  914.2000  420.0000];
nexttile(1);
plot(R_vals_analytic/2,T_R_analytic,'--','Color','k');
T_R = ones(N_rr,1)*T_cyl_0;
xrange = [0,max(R_vals/2)*1.1];
zrange = [T_air,T_cyl_0*1.1];
hold on;
plot(R_vals/2,T_R);
xline(0.5,':','R_i = 1','Color','k')
xline(ro_vx,':',sprintf('Ro = %d',ro_vx),'Color','k');
ylim(zrange);
xlim(xrange);

nexttile(2);
T_2D_plot = T_2D_init;
T_2D_plot(ind_pipe(1),ind_pipe(2)) = nan;
T_2D_plot(~cylinder_bool) = nan;
him = imagesc(T_2D_plot,"XData",[-ro_vx,ro_vx],"YData",[-ro_vx,ro_vx]);
cbar = colorbar();
clim([0,1])
cbar.Label.String = '\theta';
cbar.Label.Rotation = 0;
cbar.Label.Position = [0.5000    1.0350         0];
cbar.Label.FontWeight = 'bold';

set(him, 'AlphaData', ~isnan(T_2D_plot))
drawnow();
hold off;

if save_animation
    frame = getframe(gcf);
    writeVideo(myVideo,frame);
end
for tt = 1:N_time_steps
    Fo_cur = tt*Fo_step;
    [M_trans, RHS_trans] = transientTerm(T_2D_cell_old, time_step); % Get transient term based on old temperature field

    iii_pipe = matrix_ind_map(ind_pipe(1)+1,ind_pipe(2)+1); 
    M_trans(iii_pipe,iii_pipe) = 0; % Removing diffusion from row associated with the pipe
    RHS_trans(iii_pipe) = 0;
    M = M_trans-M_diff+M_pipe + M_BC;
    RHS = RHS_trans+RHS_BC + RHS_pipe;

    T_2D_cell_new = solvePDE(mesh_2D,M, RHS); % Solving linear equation M*T=RHS -> this basically just calls M\RHS
    % Update temperature field for FVTool transient term
    T_2D_cell_old = T_2D_cell_new;
    
    T_2D_new = T_2D_cell_new.value(2:(end-1),2:(end-1)); % updating temperatures (ignoring ghost cells)
    T_R = zeros(N_rr,1);
    counts = zeros(N_rr,1);
    for ii = 1:Nx
        x_dif = ii-0.5-center(1);
        for jj = 1:Ny
            y_dif = jj-0.5-center(2);
            rr = round(sqrt(x_dif^2+y_dif^2));
            if rr > N_rr || rr == 0
                continue
            end
            T_R(rr) = T_R(rr) + T_2D_new(ii,jj);
            counts(rr) = counts(rr)+1;
        end
    end
    T_R = T_R./counts;        

    sol = HollowCylinderConvectionTransient(Fo_cur*Ro^2,lambdas,W,E);
    T_R_analytic = T_air + sol;

    % Update plot
    if mod(tt,N_stop)==0%mod(tt/N_time_steps*100,1) == 0
        nexttile(1)
        plot(R_vals_analytic/2,T_R_analytic,'--','Color','k');
    
        xlim(xrange);ylim(zrange);
        hold on
        plot(R_vals/2,T_R);
        xline(ri_vx,':','r_i = 0.5','Color','k')
        xline(ro_vx,':',sprintf('r_o = %0.1f',ro_vx),'Color','k');
        legend({'Analytic','Finite volume method'},'Location','southeast')
        xlabel('r'); ylabel('\theta = T/T_0')
       
        hold off
        nexttile(2)
        T_2D_plot = T_2D_new;
        T_2D_plot(ind_pipe(1),ind_pipe(2)) = nan;
        T_2D_plot(~cylinder_bool) = nan;
        him = imagesc(T_2D_plot,"XData",[-ro_vx,ro_vx],"YData",[-ro_vx,ro_vx]);
        set(him, 'AlphaData', ~isnan(T_2D_plot))
        cbar = colorbar();
        clim([0,1])
        cbar.Label.String = '\theta';
        cbar.Label.FontWeight = 'bold';
        cbar.Label.Rotation = 0;
        cbar.Label.Position = [0.5000    1.0350         0];
        xlims = xlim();
        ylims = ylim();
        text(xlims(1)+0.82*(xlims(2)-xlims(1)),ylims(1)+0.05*(ylims(2)-ylims(1)),sprintf('Fo = %0.2f',Fo_cur))
        
        drawnow()
        pause(0.01);
        
        if save_animation
            frame = getframe(gcf);
            writeVideo(myVideo,frame);
        end
       
    end
end
close(myVideo)

%% END SCRIPT %%
function sol = HollowCylinderConvectionTransient(Fo,lambdas,W,E)
        sol = sum(E.*exp(-(lambdas.^2).*Fo).*W,1);
end

function W = get_W(R,lambdas,Bi_i)
    W = -(Bi_i.*Y0(lambdas) + lambdas.*Y1(lambdas)).*J0(lambdas.*R) ...
        +(lambdas.*J1(lambdas)+Bi_i*J0(lambdas)).*Y0(lambdas.*R);
end

function E = get_E(lambdas,Bi_i,Bi_o,T_1,T_f,Ro)
    intfun = @(Rint,lambda) Rint.*(T_1-T_f).*get_W(Rint,lambda,Bi_i);
    intvals = zeros(length(lambdas),1);
    for ii = 1:length(lambdas)
        intfun_cur = @(Rint) intfun(Rint,lambdas(ii));
        intvals(ii) = integral(intfun_cur,1,Ro);
    end
    E = pi^2/2*lambdas.^2.*(Bi_o*J0(lambdas*Ro)-lambdas.*J1(lambdas*Ro)).^2.*intvals...
        .*((lambdas.^2+Bi_o^2).*(Bi_i*J0(lambdas)+lambdas.*J1(lambdas)).^2 ...
        -(lambdas.^2+Bi_i^2).*(Bi_o*J0(lambdas*Ro)-lambdas.*J1(lambdas*Ro)).^2).^(-1);    
end

function val = lambda_fun(lambda,Bi_i,Bi_o,Ro)
 val = (Bi_i*J0(lambda)+lambda.*J1(lambda)).*(Bi_o*Y0(lambda*Ro)-lambda.*Y1(lambda*Ro)) ...
     - (Bi_o*J0(lambda*Ro)-lambda.*J1(lambda*Ro)).*(Bi_i*Y0(lambda)+lambda.*Y1(lambda));
end

function lambdas = get_lambdas(N,Bi_i,Bi_o,Ro)
    syms lambda
    eqn = (Bi_i*J0(lambda)+lambda*J1(lambda))*(Bi_o*Y0(lambda*Ro)-lambda*Y1(lambda*Ro))...
        - (Bi_o*J0(lambda*Ro)-lambda*J1(lambda*Ro))*(Bi_i*Y0(lambda)+lambda*Y1(lambda)) ==0;
    safety_factor = 10;
    sampling_freq = 0.1592*(Ro-1)*safety_factor; % Determined this value of 0.1592 seperately
    sampling_period = 1/sampling_freq;
    lambdas = zeros(N,1);
    n = 0;
    lambda_found = 0;
    t_prev = 1e-6; % offset it a bit
    val_prev = lambda_fun(t_prev,Bi_i,Bi_o,Ro); % First value very close to 0
    while lambda_found < N
        n = n+1;
        t_cur = t_prev+sampling_period;
        val_cur = lambda_fun(t_cur,Bi_i,Bi_o,Ro);
        if val_cur*val_prev<0
            lambda_cur = double(vpasolve(eqn,lambda,[t_prev t_cur]));
            err = lambda_fun(lambda_cur,Bi_i,Bi_o,Ro);
            if err<1e-2
                lambda_found = lambda_found +1;
                lambdas(lambda_found) = lambda_cur;
            end
        end
        val_prev = val_cur;
        t_prev = t_cur;
    end
end

%% Define the bessel functions in shorthand just for readability
function J0 = J0(x)
    J0 = besselj(0,x);
end

function J1 = J1(x)
    J1 = besselj(1,x);
end

function Y0 = Y0(x)
    Y0 = bessely(0,x);
end

function Y1 = Y1(x)
    Y1 = bessely(1,x);
end


