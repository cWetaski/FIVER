% Goal: given mass flow rate, wall temperature, length, diameter,
% compute exit temperature
clc
clear
close all

% Properties:
% For now, let's just use a circular stainless steel pipe with air:
% eps_steel = 2*10^(-3) mm https://www.engineersedge.com/fluid_flow/pipe-roughness.htm
% Pr_air = 0.71 -> functional relationship with temperature but stays within 0.7-0.74 for T in 0-1000°C at 1 bar, https://www.engineeringtoolbox.com/air-prandtl-number-viscosity-heat-capacity-thermal-conductivity-d_2009.html

%% Calculation parameters:
T_w = 60 + 273.15; % K: wall temperature
T_i = 20 + 273.15; % K: air Inlet temperature
L = 0.25; % m length of the section
N_segments = 100; % Number of segments to divide pipe into for second method
L_seg = L/N_segments;
x_edges = (0:N_segments)*L_seg;
x_centers = (x_edges(1:end-1)+x_edges(2:end))/2;
D = 0.01;
u_flow = 0.7;


%% Define properties
eps_steel = 2*10^(-3)/1000; % m: Roughness of stainless steel
p = 101325; % Pa: Pressure -> it is fair to assume pressure variation is small between inlet and oultet.
R_air = 287; % J/kg-K: Gas constant of air

air_table = readtable('air_properties.txt','Delimiter',' ');

%% Derived values    
T_b = T_i; % Initially use inlet temperature for evaluating properties
T_diff = 1e6;
num_itr = 0;
A_cs = pi/4*D.^2; % m2: Cross sectional area of the circular duct
perimiter = pi*D; % m: Perimiter
T_b = (T_i+T_w)/2; % The temperature they use for properties in the example
[mu,cp,k] = GetAirProps(T_b,air_table); % inlet properties
rho = p./R_air./T_b; % kg/m3: Density at inlet (Ideal gas law)  
m_dot = A_cs*u_flow*rho; % kg/s -> this is constant



while T_diff > 1 && num_itr < 10
    num_itr = num_itr + 1;
    T_b_prev = T_b;
    [mu,cp,k] = GetAirProps(T_b,air_table);
    rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
    u_flow = m_dot./(rho.*A_cs); % m/s: air velocity

    nu = mu./rho; % m2/s: Kinematic viscosity 
    Pr = cp.*mu./k; % Prandtl number
    Re = D*u_flow/nu; % Reynolds number; 
    f_darcy = GetDarcy(Re,eps_steel/D);
    f_fanning = f_darcy/4;
    Re_roughness = eps_steel/D*Re*sqrt(f_darcy/2); % Roughness reynolds number: Eq 7.85 from Sekulic, 2023

    if Re < 2300
        Gz_R = Re*Pr*D/L;
        Nu_inf_am = 3.657; % Constant wall temperature
        Nu_mean_lam = 3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R); % Lienhard, Eq 7.29
    end
    if Re > 2100
        Nu_inf_turb = (f_darcy/8).*(Re-1000).*Pr./(1+12.7.*(f_darcy/8).^0.5.*(Pr.^(2/3)-1)); % Eq 7.41 from lienhard/7.76 from Sekulic
        C6 = (L./D).^0.1./(Pr.^(1/6)).*(0.68+3000./Re.^0.81); % Eq 7.102 from Sekulic, 2023
        Nu_mean_turb = (1 + C6./(L./D)).*Nu_inf_turb; % Eq 7.101 from Sekulic, 2023, for thermally developing flows (note, lienhard's entrance correlation for turbulent flows is for simultaneously developing flow 
        %Nu_mean_turb = (1+0.9756/(L/D)^0.76)*Nu_inf_turb;
        %Nu_mean_turb = Nu_mean_turb*(T_b/T_w)^(0.47); % Eq 7.45 from Lienhard
        n = -log10(T_w/T_b)^(-1/4)+0.3;
        Nu_mean_turb = Nu_mean_turb*(T_w/T_b)^n; % Correction for wall temperature, from Sekulic;
        
    end

    if Re <= 2100 % Laminar, thermally developing
        Nu_mean = Nu_mean_lam;
    elseif Re >= 2300
        Nu_mean = Nu_mean_turb;
    else 
        alpha = (2300-Re)/(2300-2100); 
        Nu_mean = alpha*Nu_mean_lam + (1-alpha)*Nu_mean_turb;  
    end
   
   h_mean = k.*Nu_mean./D; % Definition of nusselt number
   T_o = T_i + (T_w-T_i).*(1 - exp(-h_mean.*perimiter.*L./(m_dot.*cp))); % Lienhard, Eq 7.57
   if num_itr == 1
       fprintf("Textbook: nu = %0.7f, k = %0.5f, rho = %0.2f,  First itr: nu = %0.7f, k = %0.5f, rho = %0.2f \n",1.69e-5,0.00270,1.13,nu,k,rho)
       fprintf("Textbook: Nu = %0.2f, h = %0.2f,  First itr:   Nu = %0.2f, h = %0.2f \n",4.27,11.5,Nu_mean,h_mean);
       
       T_o_C_first_itr = T_o-273.15; % Store first answer, as it should perfectly agree with Lienhard
   end
   Q_tot = m_dot.*cp.*(T_o-T_i);
    
   T_b = (T_o+T_i)/2;
   T_diff = abs(T_b-T_b_prev);
end
fprintf("T_b (mean, final) = %0.2f°C \n",T_b-273.15)
% Now do it by pipe segment
T_i_x = zeros(1,N_segments);
T_o_x = zeros(1,N_segments);
Nu_x = zeros(1,N_segments);
T_i_x(1) = T_i;
T_diff = 1e6;
count = 0;
T_b = T_i_x(1);
%while T_diff > 0.1 && count<100
for jj = 1:N_segments
    T_diff = 1e6;
    T_b = T_i_x(jj);
    num_itr = 0;
    while T_diff > 1 && num_itr < 10
        num_itr = num_itr+1;
        x_L = x_edges(jj);
        x_R = x_edges(jj+1);
        
        [mu,cp,k] = GetAirProps(T_b,air_table);
        %[mu,cp,k] = GetAirProps(T_b,air_table);

        rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
        %rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
        
        nu = mu./rho; % m2/s: Kinematic viscosity 
        Pr = cp.*mu./k; % Prandtl number
        
        Re_x = D*u_flow/nu; % Reynolds number;
        f_darcy = GetDarcy(Re_x,eps_steel/D); 
        f_fanning = f_darcy/4;
        if Re_x < 2300
            Gz_L = Re_x*Pr*D/x_L;
            Gz_R = Re_x*Pr*D/x_R;
            if x_L == 0
                Nu_x_lam = 3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R);
            else
                Nu_x_lam = 1/L_seg*( x_R*(3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R)) ...
                                    -x_L*(3.657/(tanh(2.264*Gz_L^(-1/3)+1.7*Gz_L^(-2/3)))+0.0499*Gz_L*tanh(1/Gz_L))); % Lienhard, Eq 7.29
            end
        end
        if Re_x > 2100
            Nu_inf_turb = (f_darcy/8).*(Re_x-1000).*Pr./(1+12.7.*(f_darcy/8).^0.5.*(Pr.^(2/3)-1)); % Eq 7.41 from lienhard/Eq 7.76 from Sekulic
            if x_L == 0
                Nu_x_turb =  Nu_inf_turb*(1 + (x_R/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81)); % Left value is undefined
            else
                Nu_x_turb = 1/L_seg*Nu_inf_turb*(  x_R*(1 + (x_R/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81)) ...
                                                  -x_L*(1 + (x_L/D)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_x^0.81))); % Eq 7.102 from sekulic, Average value from x_L to x_R       
            end
            %Nu_x_turb = Nu_x_turb*(T_b/T_w)^(0.47); % Correction for wall temperature, from Lienhard:
            n = -log10(T_w/T_b)^(-1/4)+0.3;
            Nu_x_turb = Nu_x_turb*(T_w/T_b)^n; % Correction for wall temperature, from Sekulic;
        end

        if Re_x <= 2100 % Fully laminar
            Nu_x(jj) = Nu_x_lam;
        elseif Re_x >= 2300 % validity range for Gnielinski
            Nu_x(jj) = Nu_x_turb;
        else % Take weighted combination
            disp(' Not fully turb')
            alpha = (2300-Re_x)/(2300-2100);
            Nu_x(jj) = alpha*Nu_x_lam+(1-alpha)*Nu_x_turb;
        end
        h_x = k.*Nu_x(jj)./D; % Eq 7.26 from Sekulic, 2023
        T_o_x(jj) = T_i_x(jj) + (T_w-T_i_x(jj)).*(1-exp(-h_x.*perimiter.*L_seg./(m_dot.*cp)));
        LMTD = (T_o_x(jj)-T_i_x(jj))./log((T_w-T_i_x(jj))./(T_w-T_o_x(jj))); % log mean temperature in the section
        T_b_prev = T_b;
        T_b = (T_o_x(jj)+T_i_x(jj))/2; % Estimate bulk temperature in section
        T_diff = abs(T_b-T_b_prev);
        if jj < N_segments
            T_i_x(jj+1) = T_o_x(jj);
        end
    end
end
%LMTD = (T_o_x(end)-T_i_x(1))./log((T_w-T_i_x(1))./(T_w-T_o_x(end))); % log mean temperature difference assuming constant wall temp
%T_b_from_LMTD = T_w - LMTD; 
%fprintf("%0.1f,    %0.1f \n",T_b,T_b_from_LMTD) % Can see that average along length is typically within 1 degree of LMTD
%T_diff = abs(T_b_prev-T_b);

Nu_x_mean = mean(Nu_x);
T_o_seg = T_o_x(end);
Q_tot_seg = m_dot.*cp.*(T_o_seg-T_i);

T_o_C = T_o-273.15;
T_o_seg_C = T_o_seg-273.15;

fprintf("Textbook: %0.1f°C,  Mean (first itr): %0.1f°C,  Mean: %0.1f°C,  Segmented: %0.1f°C \n",50.6,T_o_C_first_itr,T_o_C,T_o_seg_C);

