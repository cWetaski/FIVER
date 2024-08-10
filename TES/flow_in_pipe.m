% Goal: given mass flow rate, wall temperature, length, diameter,
% compute exit temperature
clc
clear
close all

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextInterpreter','latex')


% Properties:
% For now, let's just use a circular stainless steel pipe with air:
% eps_steel = 2*10^(-3) mm https://www.engineersedge.com/fluid_flow/pipe-roughness.htm
% Pr_air = 0.71 -> functional relationship with temperature but stays within 0.7-0.74 for T in 0-1000°C at 1 bar, https://www.engineeringtoolbox.com/air-prandtl-number-viscosity-heat-capacity-thermal-conductivity-d_2009.html

%% Calculation parameters:
T_w = 400 + 273.15; % K: wall temperature
T_i = 0 + 273.15; % K: air Inlet temperature
L = 1; % m length of the section
N_segments = 11; % Number of segments to divide pipe into for second method
L_seg = L/N_segments;
m_dot = 0.002; % kg/s
x_edges = (0:N_segments)*L_seg;
x_centers = (x_edges(1:end-1)+x_edges(2:end))/2;

%% Dependent variable:
D_all = linspace(0.01,0.1,100); % m;
%D_all = 0.1;
%% Define properties
roughness_steel = 2*10^(-3)/1000; % m: Roughness of stainless steel
p = 1e5; % Pa: Pressure -> it is fair to assume pressure variation is small between inlet and oultet.
R_air = 287; % J/kg-K: Gas constant of air

air_table = readtable('air_properties.txt','Delimiter',' ');

%% Derived values    
for ii = 1:length(D_all)
    T_b = T_i; % Initially use inlet temperature for evaluating properties
    T_diff = 1e6;
    num_itr = 0;
    D = D_all(ii);
    A_cs = pi/4*D.^2; % m2: Cross sectional area of the circular duct
    perimiter = pi*D; % m: Perimiter
    

    while T_diff > 1 && num_itr < 10
        num_itr = num_itr + 1;
        T_b_prev = T_b;
        [mu,cp,k] = GetAirProps(T_b,air_table);
        rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
        nu = mu./rho; % m2/s: Kinematic viscosity 
        Pr = cp.*mu./k; % Prandtl number
        
        u_flow = m_dot./(rho.*A_cs); % m/s: air velocity
     
        Re = D*u_flow/nu; % Reynolds number;
        Re_all(ii) = Re;  
        f_darcy = GetDarcy(Re,roughness_steel/D);
        f_fanning = f_darcy/4;
        Re_roughness(ii) = roughness_steel/D*Re*sqrt(f_darcy/2); % Roughness reynolds number: Eq 7.85 from Sekulic, 2023

        if Re < 2300
            Gz_R = Re*Pr*D/L;
            Nu_inf_am = 3.657; % Constant wall temperature
            Nu_mean_lam = 3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R); % Lienhard, Eq 7.29
            f_darcy_lam = f_darcy*(T_w/T_b); % from Sekulic, 2023
        end
        if Re > 2100 % Turbulent/Transition
            Nu_inf_turb = (f_darcy/8).*(Re-1000).*Pr./(1+12.7.*(f_darcy/8).^0.5.*(Pr.^(2/3)-1)); % Eq 7.41 from lienhard/7.76 from Sekulic
            C6 = (L./D).^0.1./(Pr.^(1/6)).*(0.68+3000./Re.^0.81); % Eq 7.102 from Sekulic, 2023
            Nu_mean_turb = (1 + C6./(L./D)).*Nu_inf_turb; % Eq 7.101 from Sekulic, 2023, for thermally developing flows (note, lienhard's entrance correlation for turbulent flows is for simultaneously developing flow 
            %Nu_mean_turb = (1+0.9756/(L/D)^0.76)*Nu_inf_turb;
            %Nu_mean_turb = Nu_mean_turb*(T_b/T_w)^(0.47); % Eq 7.45 from Lienhard
            n = -log10(T_w/T_b)^(-1/4)+0.3;
            Nu_mean_turb = Nu_mean_turb*(T_w/T_b)^n; % Correction for wall temperature, from Sekulic;
            f_darcy_turb = f_darcy*(T_w/T_b)^(-0.1); % from Sekulic, 2023
        end

        if Re <= 2100 % Laminar, thermally developing
            Nu_mean(ii) = Nu_mean_lam;
            f_darcy = f_darcy_lam;
        elseif Re >= 2300 % Turbulent or transition where Gnielinski correlation is valid
            Nu_mean(ii) = Nu_mean_turb;
            f_darcy = f_darcy_turb;
        else 
            alpha = (2300-Re)/(2300-2100); 
            Nu_mean(ii) = alpha*Nu_mean_lam + (1-alpha)*Nu_mean_turb;  
            f_darcy = alpha*f_darcy_lam + (1-alpha)*f_darcy_turb;
        end

        h_mean = k.*Nu_mean(ii)./D; % Definition of nusselt number
        
       T_o(ii) = T_i + (T_w-T_i).*(1 - exp(-h_mean.*perimiter.*L./(m_dot.*cp))); % Lienhard, Eq 7.57
       Q_tot(ii) = m_dot.*cp.*(T_o(ii)-T_i);
        
       T_b = (T_o(ii)+T_i)/2;
       T_diff = abs(T_b-T_b_prev);
    end
    % Now do it by pipe segment
    p_loss(ii) = L*f_darcy*rho/2*u_flow.^2/D; % Pressure drop
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
            
            u_flow = m_dot./(rho.*A_cs); % m/s: air velocity
            Re_x = D*u_flow/nu; % Reynolds number;
            Re_x_all(jj) = Re_x;
            f_darcy = GetDarcy(Re_x,roughness_steel/D); 
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
                f_darcy_lam = f_darcy*(T_w/T_b); % from Sekulic, 2023
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
                f_darcy_turb = f_darcy*(T_w/T_b)^(-0.1); % from Sekulic, 2023
            end
    
            if Re_x <= 2100 % Fully laminar
                Nu_x(jj) = Nu_x_lam;
                f_darcy = f_darcy_lam;
            elseif Re_x >= 2300 % validity range for Gnielinski
                Nu_x(jj) = Nu_x_turb;
                f_darcy = f_darcy_turb;
            else % Take weighted combination
                disp('Weighted avg')
                alpha = (2300-Re_x)/(2300-2100);
                Nu_x(jj) = alpha*Nu_x_lam+(1-alpha)*Nu_x_turb;
                f_darcy = alpha*f_darcy_lam+(1-alpha)*f_darcy_turb;
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
            p_loss_x(jj) = L_seg*f_darcy*rho/2*u_flow.^2/D; % Pressure drop in segment
        end
    end
    %LMTD = (T_o_x(end)-T_i_x(1))./log((T_w-T_i_x(1))./(T_w-T_o_x(end))); % log mean temperature difference assuming constant wall temp
    %T_b_from_LMTD = T_w - LMTD; 
    %fprintf("%0.1f,    %0.1f \n",T_b,T_b_from_LMTD) % Can see that average along length is typically within 1 degree of LMTD
    %T_diff = abs(T_b_prev-T_b);
    
    Nu_x_mean(ii) = mean(Nu_x);
    T_o_seg(ii) = T_o_x(end);
    Q_tot_seg(ii) = m_dot.*cp.*(T_o_seg(ii)-T_i);
    p_loss_seg(ii) = sum(p_loss_x);
end
f = figure;
x_locs = L_seg.*(2.*(1:N_segments)-1)/2; % center of position of each segment
            
plot(x_locs,Nu_x)
eta_pump = 0.8;
pumping_energy = p_loss.*m_dot./rho/eta_pump;



f1 = figure;
plot(D_all,Re_roughness)
yline(5,'--')
x_min = min(D_all); x_max = max(D_all); x_range = x_max-x_min;
xlim([x_min,x_max]);

y_max = max(Re_roughness);
ylim([0,y_max]);

if y_max > 5
    text(x_min + x_range*0.7,5-y_max*0.04,'Smooth regime','HorizontalAlignment','Center')
    text(x_min + x_range*0.7,5+y_max*0.04,'Transition regime','HorizontalAlignment','Center')
else
    text(x_min + x_range*0.3,y_max*0.7,'Smooth regime: (Roughness Reynolds $<$ 5)')
end
if y_max > 70
    text(x_min + x_range*0.7,70+y_max*0.04,'Rough regime','HorizontalAlignment','Center');
end

xlabel('Diameter (m)')
ylabel('Roughness Reynolds number')

f2 = figure;
hold on
yyaxis left
p1 = plot(D_all,Nu_mean,'Color',[0 0.4470 0.7410]);
p2 = plot(D_all,Nu_x_mean,'--','Color',[0 0.4470 0.7410]);
ylabel("Mean Nusselt number")


yyaxis right
plot(D_all,Re_all)
ylabel("Reynolds number")
xlabel('Diameter (m)')
legend([p1,p2],'Mean','Segments')

f3 = figure;
hold on;
xlabel('Diameter (m)')
yyaxis left
plot(D_all,T_o-273.15,'Color',[0 0.4470 0.7410]);
plot(D_all,T_o_seg-273.15,'--','Color',[0 0.4470 0.7410])
ylabel('Outlet temperature ($^{\circ}$C)')
yyaxis right
plot(D_all,Q_tot,'Color',[0.8500 0.3250 0.0980]);
plot(D_all,Q_tot_seg,'--','Color',[0.8500 0.3250 0.0980])
ylabel('Heat transfer rate (W)')

leg1 = plot(nan,nan,'-','Color','k');
leg2 = plot(nan,nan,'--','Color','k');
legend([leg1,leg2],"Mean","In segments")


f4 = figure;
hold on;
plot(D_all,p_loss,'-','Color','k'); % refer to pressure class for duct system: https://www.engineeringtoolbox.com/duct-systems-pressure-classification-d_2150.html
plot(D_all,p_loss_seg,'--','Color','k')
legend("Mean","Segment")
ylabel('Pumping pressure loss (Pa)')
xlabel('Diameter (m)')
xlim([0.05,max(D_all)])