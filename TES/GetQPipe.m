function [h_mean_x,T_b_x,Q_tot_x,T_o] = GetQPipe(T_w_x,T_air_i,m_dot,L_seg,D,roughness)
%GETQPIPE Summary of this function goes here
% All scalar variables except T_w, which should be a 1D vector.
%   Detailed explanation goes here
    
    p = 101325; % Pa;
    air_table = readtable('air_properties.txt','Delimiter',' '); % Table of air properties vs T at p=1atm
    [~,cp_i,~] = GetAirProps(T_air_i,air_table); % Get initial cp

    R_air = 287; % J/kg-K: Gas constant of air
    
    Nx = length(T_w_x);
    x_edges = (0:Nx)*L_seg;
    A_cs = pi/4*D.^2; % m2: Cross sectional area of the circular duct
    perimiter = pi*D; % m: Perimiter
    
    T_i_x = zeros(1,Nx);
    T_o_x = zeros(1,Nx);
    Nu_x = zeros(1,Nx);
    T_b_x = zeros(1,Nx);
    T_i_x(1) = T_air_i;
    %while T_diff > 0.1 && count<100
    Delta_T = 0;
    for ii = 1:Nx
        T_diff = 1e6;
        T_b = T_i_x(ii)+Delta_T/2; % Approximate first guess using Delta T from previous segment
        num_itr = 0;
        T_w = T_w_x(ii);
        while T_diff > 1 && num_itr < 10
            num_itr = num_itr+1;
            x_L = x_edges(ii);
            x_R = x_edges(ii+1);
            try
            [mu,cp,k] = GetAirProps(T_b,air_table);
            %[mu,cp,k] = GetAirProps(T_b,air_table);
            catch e
                debug = 0;
            end
            rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
            %rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
            
            nu = mu./rho; % m2/s: Kinematic viscosity 
            Pr = cp.*mu./k; % Prandtl number
            
            u_flow = m_dot./(rho.*A_cs); % m/s: air velocity
            Re_x = D*u_flow/nu; % Reynolds number;
            f_darcy = GetDarcy(Re_x,roughness/D); 
            if Re_x < 2300
                Gz_L = Re_x*Pr*D/x_L;
                Gz_R = Re_x*Pr*D/x_R;
                if x_L == 0
                    Nu_x_lam = 3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R);
                else
                    Nu_x_lam = 1/L_seg*( x_R*(3.657/(tanh(2.264*Gz_R^(-1/3)+1.7*Gz_R^(-2/3)))+0.0499*Gz_R*tanh(1/Gz_R)) ...
                                        -x_L*(3.657/(tanh(2.264*Gz_L^(-1/3)+1.7*Gz_L^(-2/3)))+0.0499*Gz_L*tanh(1/Gz_L))); % Lienhard, Eq 7.29
                end
                if T_w > T_b
                    f_darcy_lam = f_darcy*(T_w/T_b); % from Sekulic, 2023
                else
                    f_darcy_lam = f_darcy*(T_w/T_b)^0.81; % Sekulic, 2023
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
                if T_w > T_b
                    n = -log10(T_w/T_b)^(-1/4)+0.3;
                    Nu_x_turb = Nu_x_turb*(T_w/T_b)^n; % Correction for wall temperature, from Sekulic, only for heating
                end
                f_darcy_turb = f_darcy*(T_w/T_b)^(-0.1); % from Sekulic, 2023 (for both heating and cooling)
            end
            if Re_x <= 2100 % Fully laminar
                Nu_x(ii) = Nu_x_lam;
                f_darcy = f_darcy_lam;
            elseif Re_x >= 2300 % validity range for Gnielinski
                Nu_x(ii) = Nu_x_turb;
                f_darcy = f_darcy_turb;
            else % Take weighted combination
                disp('Warning: Transition')
                alpha = (2300-Re_x)/(2300-2100);
                Nu_x(ii) = alpha*Nu_x_lam+(1-alpha)*Nu_x_turb;
                f_darcy = alpha*f_darcy_lam+(1-alpha)*f_darcy_turb;
            end
            h_mean = k.*Nu_x(ii)./D; % Eq 7.26 from Sekulic, 2023
            T_o_x(ii) = T_i_x(ii) + (T_w-T_i_x(ii)).*(1-exp(-h_mean.*perimiter.*L_seg./(m_dot.*cp)));
            T_b_prev = T_b;
            T_b = (T_o_x(ii)+T_i_x(ii))/2; % Estimate bulk temperature in section
            T_diff = abs(T_b-T_b_prev);
            if ii < Nx
                T_i_x(ii+1) = T_o_x(ii);
            end
            p_loss_x(ii) = L_seg*f_darcy*rho/2*u_flow.^2/D; % Pressure drop in segment
            h_mean_x(ii) = h_mean;
        end
        T_b_x(ii) = T_b;
        Delta_T = T_o_x(ii)-T_i_x(ii);
        Q_tot_x(ii) = m_dot.*cp.*(T_o_x(ii)-T_i_x(ii)); % Heat transfer rate in the section;
    end
    %LMTD = (T_o_x(end)-T_i_x(1))./log((T_w-T_i_x(1))./(T_w-T_o_x(end))); % log mean temperature difference assuming constant wall temp
    %T_b_from_LMTD = T_w - LMTD; 
    %fprintf("%0.1f,    %0.1f \n",T_b,T_b_from_LMTD) % Can see that average along length is typically within 1 degree of LMTD
    %T_diff = abs(T_b_prev-T_b);
    
    Nu_x = mean(Nu_x);
    T_o = T_o_x(end);
    Q_tot = m_dot.*(T_o_x(end)*cp-T_air_i*cp_i);
    p_loss = sum(p_loss_x);
end

