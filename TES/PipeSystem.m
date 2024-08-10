classdef (Abstract) PipeSystem
    %PIPESYSTEM Summary of this class goes here
    % Constructors and abstract methods getPipeGeometry, getPipeEquation, and getWallTemperature, get_f_Darcy are implemented in each subclass (CirclePipe and RectanglePipe)
    
    properties
        x_loc;  % x location of center (voxel coords)
        y_loc;  % y location of center (voxel coords)
        D_hydraulic; % Hydraulic diameter (m)
        A_cs; % Cross sectional area (m2)
        VS_pipe_2D; % 2D Matrix of the entire pipe.
        VS_boundaries_2D; %
        inds_pipe_boundary; % 2 column vectors [X,Y]: Pairs correspond to the position in the voxel space (not considering ghost cells) of the boundary elements of the pipe 
        inds_pipe_interior; % 2 column vectors [X,Y]: Pairs correspond to the position in the voxel space (not considering ghost cells) of the interior pipe elements
        boundary_perimiter_span; % 1D Vector, Each element contains the perimiter (in voxel units) associated with the corresponding element of inds_pipe_boundary
        size_VS; %
        vx_scale; % (1x3 double) m/vx (everything assumes vx_scale(1)==vx_scale(2)!!)
        roughness; % m material property
        fluid_property_table; % Table
        S_wall_z; % Source term at the wall along the axis
        last_T_wall_z; % T_wall_z stored from last time it was computed
        last_T_fluid_z; % T_fluid_z stored from last time it was computed
        last_T_fluid_o; % 
        last_m_dot; % Last inlet Reynolds number
    end
    
    methods    
        
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
        %% getHeatTransferCoefficient
        function [h_mean_z,T_fluid_z,Q_flow,T_fluid_o,p_loss,Re_in] = getHeatTransferCoefficient(obj,T_w_z,T_fluid_i,m_dot)
            %GETQPIPE Summary of this function goes here
            % All scalar variables except T_w_z, which should be a 1D vector.
            %   T_w_z (1D vector [K]): Wall temperature at each z-location
            %   T_fluid_i   (scalar [K]): Inlet fluid temperature
            %   m_dot (scalar [kg/s]: Mass flow rate
            %   Re_in (Reynolds number of input)
            if m_dot == 0
                T_fluid_z = T_w_z;
                h_mean_z = zeros(size(T_w_z));
                Q_flow = 0;
                T_fluid_o = T_fluid_i;
                p_loss = 0;
                return
            end
            D_pipe = obj.D_hydraulic;
            Re_min = 1e6; Re_max = 0;
            obj.last_m_dot = m_dot;
            p = 101325; % Pa; nominal pressure
            
            L_seg = obj.vx_scale(3); % Length of each "segment" is always 1 voxel
            R_air = 287; % J/kg-K: Gas constant of air

            Nz = length(T_w_z);
            z_edges = (0:Nz)'*L_seg;
            perimiter = 4*obj.A_cs/D_pipe; % m: Perimiter

            [mu,cp_b] = obj.getFluidProperties(T_fluid_i); % Initial cp is used as the first guess for cp_b in the first section
            rho = p./R_air./T_fluid_i; % kg/m3: Density (Ideal gas law)
            %rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
            
            nu = mu./rho; % m2/s: Kinematic viscosity 
            u_flow = m_dot/(rho*obj.A_cs);
            Re_in = u_flow*obj.D_hydraulic/nu;
            
            T_i_z = zeros(Nz,1);
            T_o_z = zeros(Nz,1);
            Nu_z = zeros(Nz,1);
            T_fluid_z = zeros(Nz,1);
            p_loss_z = zeros(Nz,1);
            h_mean_z = zeros(Nz,1);
            Q_tot_z = zeros(Nz,1);

            T_i_z(1) = T_fluid_i;
            %while T_diff > 0.1 && count<100
            Delta_T = 0;
            for ii = 1:Nz
                T_diff = 1e6;
                T_b = T_i_z(ii)+Delta_T/2; % Approximate first guess using Delta T from previous segment
                num_itr = 0;
                T_w = T_w_z(ii);
                while T_diff > 0.1 && num_itr < 100
                    num_itr = num_itr+1;
                    z_lower = z_edges(ii);
                    z_upper = z_edges(ii+1);
                    
                    [mu,~,k] = obj.getFluidProperties(T_b); % cp_b is set at the end of the loop, in the first iteration, the cp_b from the previous segment is used
                    rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
                    %rho = p./R_air./T_b; % kg/m3: Density (Ideal gas law)
                    
                    nu = mu./rho; % m2/s: Kinematic viscosity 
                    Pr = cp_b.*mu./k; % Prandtl number
                    
                    u_flow = m_dot./(rho.*obj.A_cs); % m/s: air velocity
                    Re_x = D_pipe*u_flow/nu; % Reynolds number;
                    if Re_x < Re_min
                        Re_min = Re_x;
                    end
                    if Re_x > Re_max
                        Re_max = Re_x;
                    end
                    
                    if Re_x < 4000
                        if Re_x <= 2300
                            Re_lam = Re_x;
                        else
                            Re_lam = 2300;
                        end
                        z_star_L = z_lower/Re_lam*Pr*D_pipe;
                        z_star_U = z_upper/Re_lam*Pr*D_pipe;
                        Nu_m_seg_lam = obj.getNu_m_lam(z_star_L,z_star_U);
                        f_Darcy_lam = obj.get_f_Darcy(Re_lam,obj.roughness)*(T_w/T_b); % Mills, 2017
                    end
                    if Re_x > 2300
                        
                        if Re_x >= 4000
                            Re_turb = Re_x;
                        else
                            Re_turb = 4000;
                        end
                        f_Darcy_sm = obj.get_f_Darcy(Re_turb,0);
                        f_Darcy_rough = obj.get_f_Darcy(Re_turb,obj.roughness);
                        Nu_inf_turb_sm = (f_Darcy_sm/8).*(Re_turb-1000).*Pr./(1+12.7.*(f_Darcy_sm/8).^0.5.*(Pr.^(2/3)-1)); % Eq 7.41 from lienhard/Eq 7.76 from Sekulic
                        if f_Darcy_rough/f_Darcy_sm < 3 % Correction for roughness (MIT-FIRES)
                            Nu_inf_turb = Nu_inf_turb_sm*(f_Darcy_rough/f_Darcy_sm)^(0.68*Pr^0.215);
                        else
                            Nu_inf_turb = Nu_inf_turb_sm*2;
                        end
                        if z_lower == 0
                            Nu_m_seg_turb =  Nu_inf_turb*(1 + (z_upper/D_pipe)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_turb^0.81));
                        else
                            Nu_m_seg_turb = 1/L_seg*Nu_inf_turb*(  z_upper*(1 + (z_upper/D_pipe)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_turb^0.81)) ...
                                                              -z_lower*(1 + (z_lower/D_pipe)^(-0.9)/Pr^(1/6)*(0.68+3000/Re_turb^0.81))); % Eq 7.102 from sekulic, Average value from x_L to x_R       
                        end
                        
                        if T_w > T_b
                            n = -0.55; % Mills 2017
                            m = -0.2;
                        else
                            n = 0; % Mills 2017
                            m = -0.1;
                        end
                        Nu_m_seg_turb = Nu_m_seg_turb*(T_w/T_b)^n; % Correction for wall temperature, from Pethukov, 1970
                        f_Darcy_turb = f_Darcy_rough*(T_w/T_b)^m; % from Mills 2017
                    end
                    if Re_x <= 2300 % Fully laminar
                        Nu_z(ii) = Nu_m_seg_lam;
                        f_Darcy = f_Darcy_lam;
                    elseif Re_x >= 4000 % Fully turbulent
                        Nu_z(ii) = Nu_m_seg_turb;
                        f_Darcy = f_Darcy_turb;
                    else % Take weighted combination
                        %disp('Warning: Transition')
                        alpha = (Re_x-2300)/(4000-2300); % Take simple weighting between the two functions
                        
                        Nu_z(ii) = (1-alpha)*Nu_m_seg_lam+ alpha*Nu_m_seg_turb;
                        f_Darcy = (1-alpha)*f_Darcy_lam+ alpha*f_Darcy_turb;
                    end
                    h_mean = k.*Nu_z(ii)./D_pipe; % Eq 7.26 from Sekulic, 2023
                    T_o_z(ii) = T_i_z(ii) + (T_w-T_i_z(ii)).*(1-exp(-h_mean.*perimiter.*L_seg./(m_dot.*cp_b)));
                   
                    T_b_prev = T_b;
                    T_b = (T_o_z(ii)-T_i_z(ii))/log(T_o_z(ii)/T_i_z(ii)); % Log mean temperature
                    
                    T_diff = abs(T_b-T_b_prev);
                    if ii < Nz
                        T_i_z(ii+1) = T_o_z(ii);
                    end
                    
                end
                T_fluid_z(ii) = T_w-m_dot/h_mean/perimiter/L_seg*cp_b*(T_o_z(ii)-T_i_z(ii));
                Delta_T = T_o_z(ii)-T_i_z(ii);
                p_loss_z(ii) = L_seg*f_Darcy*rho/2*u_flow.^2/D_pipe; % Pressure drop in segment (Pa)
                h_mean_z(ii) = h_mean;
                Re_x_all(ii) = Re_x;
                cp_b = integral(@(T) obj.getFluidcp(T),T_i_z(ii),T_o_z(ii))/Delta_T; % exact average specific heat across the temperature difference
                Q_tot_z(ii) = m_dot.*cp_b*Delta_T; % (W) Heat transfer rate in the section ;
                k_z(ii) = k;
                
            end
            %test = sum(h_mean_z.*(T_w_z-T_fluid_z)*perimiter*L_seg);
            %fprintf('Q_air: %0.4f kW, Q_pipe: %0.4f kW \n',sum(Q_tot_z)/1000,test/1000)
            %LMTD = (T_o_x(end)-T_i_x(1))./log((T_w-T_i_x(1))./(T_w-T_o_x(end))); % log mean temperature difference assuming constant wall temp
            %T_b_from_LMTD = T_w - LMTD; 
            %fprintf("%0.1f,    %0.1f \n",T_b,T_b_from_LMTD) % Can see that average along length is typically within 1 degree of LMTD
            %T_diff = abs(T_b_prev-T_b);
            
            %Nu_mean = mean(Nu_z);
            T_fluid_o = T_o_z(end);
            C_int = integral(@(T) obj.getFluidcp(T),T_fluid_i,T_fluid_o);    
            Q_flow2 = m_dot*C_int;
            Q_flow = sum(Q_tot_z);
            fprintf("Q1 = %0.4f kW, Q2 = %0.4f kW\n",Q_flow/1000,Q_flow2/1000)
            %Q_tot = m_dot.*(T_o_x(end)*cp-T_air_i*cp_i);
            p_loss = sum(p_loss_z);
            %fprintf("Re_{min}: %0.0f, Re_{max}: %0.0f\n",Re_min,Re_max);
        end
        function [mu,cp,k,rho] = getFluidProperties(obj,T)
            % Interpolate values from appropriatetly defined table (table values in C)
            mu = interp1(obj.fluid_property_table.T,obj.fluid_property_table.mu,T-273.15)*1e-5; % kg/m-s
            cp = interp1(obj.fluid_property_table.T,obj.fluid_property_table.cp,T-273.15); % J/kg-K
            k = interp1(obj.fluid_property_table.T,obj.fluid_property_table.k,T-273.15); % W/m-K
            rho = interp1(obj.fluid_property_table.T,obj.fluid_property_table.rho,T-273.15);
        end
        function cp = getFluidcp(obj,T)
            [~,cp] = obj.getFluidProperties(T);
        end
    end
    methods (Abstract)
        getPipeGeometry(obj), 
        getPipeEquation(obj),
        getWallTemperature(obj),
        get_f_Darcy(obj)
    end
    methods (Static,Abstract)
        getNu_x_lam(x_star)
        getNu_m_lam(x_star)
    end

end

