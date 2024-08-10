clear all
close
clc

a = CirclePipe(10,10,2,0.003/1000,[1/100,1/100,1],[20,20,1],readtable("air_properties.txt"));
T_w = 100+273.15;
T_fluid_i = 273.15;

N_m = 10;
m_dot = 10.^(linspace(-3.5,-2.5,N_m));

%for ii = 1:N_m
%[h_mean(ii),~,~,~,~,Re_in(ii)] = a.getHeatTransferCoefficient(T_w,T_fluid_i,m_dot(ii));
%end

%loglog(Re_in,h_mean(ii))

Re_lam = 1000;
Pr = 0.7;
D_h = 1/100*2;

x = linspace(0,1,N_m+1);
x(1) = [];

for i = 1:N_m
    if i == 1
        xl = 0;
    else
        xl = x(i-1);
        
    end
    xr = x(i);
    Gz_l(i) = Re_lam*Pr*D_h./xl;
    Gz_r = Re_lam*Pr*D_h./xr;
    
    Nu_x(i) = getNu_x_lam(Re_lam*Pr*D_h/x(i));
    Nu_m(i) = 3.657/(tanh(2.264*Gz_r^(-1/3)+1.7*Gz_r^(-2/3)))+0.0499*Gz_r*tanh(1/Gz_r);
    if Gz_r <= 1000 && Gz_l(i) > 1000
        x1k = Re_lam*Pr*D_h/1000;
        Nu_x2(i) = (integral(@(x) getNu_x_lam(Re_lam*Pr*D_h./x),xl,x1k) ...
             + integral(@(x) getNu_x_lam(Re_lam*Pr*D_h./x),x1k,xr))/(xr-xl);
    else
        Nu_x2(i) = integral(@(x) getNu_x_lam(Re_lam*Pr*D_h./x),xl,xr)/(xr-xl);
    end
end
loglog(x,Nu_x)
hold on
x_av = ([0,x(1:end-1)]+x)/2;
loglog(x,Nu_x2)
loglog(x,Nu_m)

integral(@(x) getNu_x_lam(Re_lam*Pr*D_h./x),0,1,'ArrayValued',false)



function Nu_x = getNu_x_lam(Gz)
    % x_star is inverse Graetz number
    % Lienhard, Eq 7.28
    if Gz <= 1000
        Nu_x = 3.657+0.2362.*Gz.^0.488.*exp(-57.2./Gz);
    else
        Nu_x = 1.077*Gz.^(1/3)-0.7;
    end
end


% 
% 
% T_prev = [1000+273.15;1000+273.15];
% 
% k = 1;
% dx = 0.1;
% dt = 0.01;
% h = 10;
% T_inf = 50+273.15;
% rho = 1;
% cp = 1;
% 
% 
% mesh = createMesh1D(2,2*dx);
% BC = createBC(mesh);
% 
% alpha = createCellVariable(mesh,k/(rho*cp));
% alpha_face = harmonicMean(alpha);
% M_diff = diffusionTerm(alpha_face);
% [M_BC,RHS_BC] = boundaryCondition(BC);
% M_BC(4,:) = 0;
% 
% M_BC(4,4) = -k/dx-h/2;
% M_BC(4,3) = k/dx-h/2;
% RHS_BC(4) = -h*T_inf;
% T_old = createCellVariable(mesh,T_prev);
% [M_trans,RHS_trans] = transientTerm(T_old,dt);
% 
% M = M_BC-M_diff+M_trans;
% RHS = RHS_BC+RHS_trans;
% 
% E_prev = T_prev*dx*rho*cp;
% T_next = M\RHS;
% E_next = T_next(2:end-1)*dx*rho*cp;
% E_change = sum(E_next-E_prev);
% Q_cfd = -E_change/dt
% Q_actual_1 = h*(T_prev(end)-T_inf)
% T_w_2 = mean(T_next(end-1:end));
% T_w_mean = mean([T_w_2,T_prev(end)]);
% Q_actual_2 = h*(T_w_2-T_inf)
% Q_actual_3 = h*(T_w_mean-T_inf)
% 
% vx_scale = [1,1,1];
% roughness = 3/1000;
% [h_mean_z] = getHeatTransferCoefficient(T_w_2,T_inf,1,0.5,vx_scale,roughness)
% 
% 
% 
% 
% function [h_mean_z,T_fluid_z,Q_flow,T_fluid_o,p_loss] = getHeatTransferCoefficient(T_w_z,T_fluid_i,m_dot,radius,vx_scale,roughness)
%             GETQPIPE Summary of this function goes here
%             All scalar variables except T_w_z, which should be a 1D vector.
%               T_w_z (1D vector [K]): Wall temperature at each z-location
%               T_fluid_i   (scalar [K]): Inlet fluid temperature
%               m_dot (scalar [kg/s]: Mass flow rate
%               Re_in (Reynolds number of input)
%             if m_dot == 0
%                 T_fluid_z = T_w_z;
%                 h_mean_z = zeros(size(T_w_z));
%                 Q_flow = 0;
%                 T_fluid_o = T_fluid_i;
%                 p_loss = 0;
%                 return
%             end
% 
%             Re_min = 1e6; Re_max = 0;
%             p = 101325; % Pa; nominal pressure
% 
%             D_pipe = radius*vx_scale(1)*2;
%             L_seg = vx_scale(3); % Length of each "segment" is always 1 voxel
%             R_air = 287; % J/kg-K: Gas constant of air
% 
%             Nz = length(T_w_z);
%             z_edges = (0:Nz)'*L_seg;
%             A_cs = pi/4*D_pipe.^2; % m2: Cross sectional area of the circular duct
%             perimiter = pi*D_pipe; % m: Perimiter
% 
%             [mu,cp,k] = getFluidProperties(T_fluid_i,'air_properties.txt');
%             rho = p./R_air./T_fluid_i; % kg/m3: Density (Ideal gas law)
%             rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
% 
%             nu = mu./rho; % m2/s: Kinematic viscosity 
%             u_flow = m_dot/(rho*A_cs);
%             Re_in = u_flow*D_pipe/nu;
% 
%             T_i_x = zeros(Nz,1);
%             T_o_x = zeros(Nz,1);
%             Nu_z = zeros(Nz,1);
%             T_fluid_z = zeros(Nz,1);
%             p_loss_z = zeros(Nz,1);
%             h_mean_z = zeros(Nz,1);
%             Q_tot_z = zeros(Nz,1);
% 
%             T_i_x(1) = T_fluid_i;
%             while T_diff > 0.1 && count<100
%             Delta_T = 0;
%             for ii = 1:Nz
%                 T_diff = 1e6;
%                 T_b = T_i_x(ii)+Delta_T/2; % Approximate first guess using Delta T from previous segment
%                 num_itr = 0;
%                 T_w = T_w_z(ii);
%                 while T_diff > 1 && num_itr < 10
%                     num_itr = num_itr+1;
%                     z_lower = z_edges(ii);
%                     z_upper = z_edges(ii+1);
% 
%                     [mu,cp,k] = getFluidProperties(T_b,'air_properties.txt');
%                     rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
%                     rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
% 
%                     nu = mu./rho; % m2/s: Kinematic viscosity 
%                     Pr = cp.*mu./k; % Prandtl number
% 
%                     u_flow = m_dot./(rho.*A_cs); % m/s: air velocity
%                     Re_x = D_pipe*u_flow/nu; % Reynolds number;
%                     if Re_x < Re_min
%                         Re_min = Re_x;
%                     end
%                     if Re_x > Re_max
%                         Re_max = Re_x;
%                     end
%                     f_darcy = get_fDarcy(Re_x,roughness,radius,vx_scale);
% 
%                     if Re_x < 4000
%                         Gz_L = Re_x*Pr*D_pipe/z_lower;
%                         Gz_R = Re_x*Pr*D_pipe/z_upper;
%                         if z_lower == 0
%                             Nu_x_lam = 3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R);
%                         else
%                             Nu_x_lam = 1/L_seg*( z_upper*(3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R)) ...
%                                                 -z_lower*(3.657/(tanh(2.264*Gz_L^(-1/3)+1.7*Gz_L^(-2/3)))+0.0499*Gz_L*tanh(1/Gz_L))); % Lienhard, Eq 7.29
%                         end
%                         if T_w > T_b
%                             f_darcy_lam = f_darcy*(T_w/T_b); % from Sekulic, 2023
%                         else
%                             f_darcy_lam = f_darcy*(T_w/T_b)^0.81; % Sekulic, 2023
%                         end
%                     end
%                     if Re_x > 2100
%                         f_darcy_sm = get_fDarcy(Re_x,0,radius,vx_scale);
%                         Nu_inf_turb_sm = (f_darcy_sm/8).*(Re_x-1000).*Pr./(1+12.7.*(f_darcy_sm/8).^0.5.*(Pr.^(2/3)-1)); % Eq 7.41 from lienhard/Eq 7.76 from Sekulic
%                         if f_darcy/f_darcy_sm < 3 % Correction for roughness
%                             Nu_inf_turb = Nu_inf_turb_sm*(f_darcy/f_darcy_sm)^(0.68*Pr^0.215);
%                         else
%                             Nu_inf_turb = Nu_inf_turb_sm*2;
%                         end
%                         if z_lower == 0
%                             Nu_x_turb =  Nu_inf_turb*(1 + (z_upper/D_pipe)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81)); % Left value is undefined
%                         else
%                             Nu_x_turb = 1/L_seg*Nu_inf_turb*(  z_upper*(1 + (z_upper/D_pipe)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81)) ...
%                                                               -z_lower*(1 + (z_lower/D_pipe)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81))); % Eq 7.102 from sekulic, Average value from x_L to x_R       
%                         end
%                         if T_w > T_b
%                             n = -log10(T_w/T_b)^(-1/4)+0.3;
%                             Nu_x_turb = Nu_x_turb*(T_w/T_b)^n; % Correction for wall temperature, from Sekulic, only for heating
%                         end
%                         f_darcy_turb = f_darcy*(T_w/T_b)^(-0.1); % from Sekulic, 2023 (for both heating and cooling)
%                     end
%                     if Re_x <= 2100 % Fully laminar
%                         Nu_z(ii) = Nu_x_lam;
%                         f_darcy = f_darcy_lam;
%                     elseif Re_x >= 4000 % Fully turbulent
%                         Nu_z(ii) = Nu_x_turb;
%                         f_darcy = f_darcy_turb;
%                     else % Take weighted combination
%                         disp('Warning: Transition')
%                         alpha = ((4000-Re_x)/(4000-2100)); % Take simple weighting between the two functions
% 
%                         Nu_z(ii) = alpha*Nu_x_lam+(1-alpha)*Nu_x_turb;
%                         f_darcy = alpha*f_darcy_lam+(1-alpha)*f_darcy_turb;
%                     end
%                     h_mean = k.*Nu_z(ii)./D_pipe; % Eq 7.26 from Sekulic, 2023
%                     T_o_x(ii) = T_i_x(ii) + (T_w-T_i_x(ii)).*(1-exp(-h_mean.*perimiter.*L_seg./(m_dot.*cp)));
%                     T_b_prev = T_b;
%                     T_b = (T_o_x(ii)+T_i_x(ii))/2; % Estimate bulk temperature in section
%                     T_diff = abs(T_b-T_b_prev);
%                     if ii < Nz
%                         T_i_x(ii+1) = T_o_x(ii);
%                     end
%                     p_loss_z(ii) = L_seg*f_darcy*rho/2*u_flow.^2/D_pipe; % Pressure drop in segment (Pa)
%                     h_mean_z(ii) = h_mean;
%                 end
%                 T_fluid_z(ii) = T_b;
%                 Delta_T = T_o_x(ii)-T_i_x(ii);
%                 Q_tot_z(ii) = m_dot.*cp.*(T_o_x(ii)-T_i_x(ii)); % Heat transfer rate in the section;
%             end
%             LMTD = (T_o_x(end)-T_i_x(1))./log((T_w-T_i_x(1))./(T_w-T_o_x(end))); % log mean temperature difference assuming constant wall temp
%             T_b_from_LMTD = T_w - LMTD; 
%             fprintf("%0.1f,    %0.1f \n",T_b,T_b_from_LMTD) % Can see that average along length is typically within 1 degree of LMTD
%             T_diff = abs(T_b_prev-T_b);
% 
%             Nu_mean = mean(Nu_z);
%             T_fluid_o = T_o_x(end);
%             Q_flow = sum(Q_tot_z);
%             Q_tot = m_dot.*(T_o_x(end)*cp-T_air_i*cp_i);
%             p_loss = sum(p_loss_z);
%             fprintf("Re_{min}: %0.0f, Re_{max}: %0.0f\n",Re_min,Re_max);
%         end
%         function [mu,cp,k,rho] = getFluidProperties(T,path)
%             fluid_property_table = readtable(path);
%             Interpolate values from appropriatetly defined table (table values in C)
%             mu = interp1(fluid_property_table.T,fluid_property_table.mu,T-273.15)*1e-5; % kg/m-s
%             cp = interp1(fluid_property_table.T,fluid_property_table.cp,T-273.15); % J/kg-K
%             k = interp1(fluid_property_table.T,fluid_property_table.k,T-273.15); % W/m-K
%             rho = interp1(fluid_property_table.T,fluid_property_table.rho,T-273.15);
%         end
%         function f_Darcy = get_fDarcy(Re,roughness,radius,vx_scale)
%             % Compute darcy friction factor
%             Re: Reynolds number
%             using the weighting method given by Avci (2019), but using the Haaland correlation for the turbulent friction factor
%             rel_roughness = roughness/(radius*vx_scale(1)); % Relative roughness = epsilon/D
% 
%             f_t = (-1.8*log10((rel_roughness/3.7).^1.11+6.9./Re)).^(-2); % Haaland correlation for turbulent region
%             f_lam = 64/Re; 
% 
%             alpha = exp(-(Cm*Re/2560)^8); Method of Avci
% 
%             f = alpha*f_lam + (1-alpha)*f_t;
%             if Re < 2100
%                 f_Darcy = f_lam;
%             elseif Re > 4000
%                 f_Darcy = f_t;
%             else
%                 alpha = (4000-Re)/(4000-2100); % Simple linear weighting of f_lam and f_t;
%                 f_Darcy = f_lam*alpha + (1-alpha)*f_t;
%             end
%         end
% 
% 
