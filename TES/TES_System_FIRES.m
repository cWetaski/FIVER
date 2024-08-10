clear
close all
clc

%% USER INPUT PARAMETERS
% Geometry
N_grid_XY = 40; % Number of gridcells across the diameter
N_grid_Z = 2;
xy_length_TES = 0.2; % [m]: x-direction
z_height_TES = 18; % [m]: z-direction

% Discharge pipe parameters
% Number of pipes in the each direction in the domain (pipe axis always in z-direction)
% pipes are automatically generated, evenly spaced based on N_pipe_side
% If the automatically generated locations results in a pipe and heater at the same location, the pipe is not generated
N_pipe_side = 1;
pipe_shape = "square"; % "circle","rectangle"
r_pipe = 0.25*xy_length_TES;
xy_width_pipe = 0.5*xy_length_TES;
pipe_roughness = 3/1000; % m: 0.0015mm for glass, 3mm for alumina source: https://www.engineersedge.com/fluid_flow/pipe-roughness.htm
T_air_i = 20 + 273.15; % Kelvin

% TES Parameters - borosilicate glass SHOTT-BK7 - This is now handled by TES_material class!
TES_material_name = "Alumina"; % SHOTT_BK7_LEE, JGS2, Alumina

ns = 2; % Neighborhood size

% Charging parameters 
% Resistive heating wire https://www.kanthal.com/en/products/material-datasheets/wire/resistance-heating-wire-and-resistance-wire/kanthal-a-1/
N_heater_side = 2; % N_charge_side cannot equal N_pipe side
T_max_heater = 1400+273.15; % [K]: max temperature of charging surface
r_heater = 4/1000/2; % [m]: 4 mm diameter wire
r_void_heater = r_heater*2; % larger void reduces surface loading and therefore maximum temperature during heating (subject to discretization error)
L_frac_heater = 1; % charging rod fraction of TES length (assume it is centered in the z_direction TES)
emissivity_heater = 0.7; % actual emissivity is 0.7

% Simulation parameters
T_TES_min = 100 + 273.15;
T_TES_i = 1200 + 273.15; % Kelvin: Initial temperature of the TES
T_TES_max = T_TES_i;
N_time_steps = [200,80,180,120,120,180,180];
time_step_size = [5,10,10,15,15,20,20]; % seconds
N_rays = [10e6,5e6,4e6,4e6,4e6,3e6,2e6];
% For alumina (no participating media radiation allows bigger time steps and fewer rays)
 N_time_steps = [180,360,72];
time_step_size = [5,25,50];
 N_rays = 1e6;
end_times = cumsum(N_time_steps.*time_step_size)
heater_power_fraction = [0];
nominal_discharge_power_fraction = [1/5];%[1/5,1/5,1/5,1/5,1/5,1/5,1/5,1/5,1/5]; % 1/h: Fraction of storage capacity

T_process = 500 + 273.15; % Kelvin
z_sample = N_grid_Z;
N_par_workers = 8;
display_bool = true;
N_module = 323;

other_notes = "discharge_k5.7"
if pipe_shape == "circle"
    pipe_str = sprintf("circle_Dhpipe%dmm",round(2*r_pipe*1000))
else
    pipe_str = sprintf("square_Dhpipe%dmm",round(xy_width_pipe*1000))
end
save_filename = sprintf("%0.0f_minutes_%s_xy%dmm_z%dm_%s_domain%d_%d_%d_%s",end_times(end)/60,TES_material_name,xy_length_TES*1000,z_height_TES,pipe_str,N_grid_XY,N_grid_XY,N_grid_Z,other_notes)

%%%%% END USER DEFINED PARAMETERS %%%%%
%% Check inputs
if N_heater_side == N_pipe_side
    error("N_heater_" + ...
        "side cannot equal N_pipe_side")
end
total_time = sum(N_time_steps.*time_step_size)

if ~isempty(gcp('nocreate'))
    delete(gcp)
end
parpool('Threads',N_par_workers)

%% Derived parameters
vx_scale(1:2) = xy_length_TES/N_grid_XY; % side length of each grid cell
vx_scale(3) = z_height_TES/N_grid_Z;
r_vx_pipe = r_pipe/vx_scale(1);
xy_width_vx_pipe = xy_width_pipe/(vx_scale(1)); % Code assumes vx_scale(1)==vx_scale(2);
%r_vx_pipe = 0.5;
%r_pipe = r_vx_pipe*vx_scale(1);
D_pipe = r_pipe*2;

void_fraction = (pi*r_pipe^2)/(xy_length_TES*xy_length_TES);
volumetric_fraction = 1-void_fraction

size_VS = [N_grid_XY,N_grid_XY,N_grid_Z];
VS_T_init = zeros(size_VS(1),size_VS(2),size_VS(3));



switch TES_material_name
    case "JGS2"
        spectral_band_edges = [0.22,0.27,2,2.5,3,3.5,4.5,1000]; % good band edges for JGS2
    case "SHOTT_BK7_LEE"
        spectral_band_edges = [0.26,0.38,1.6,2.03,2.68,1000]; % Good band edges for SHOTT BK7 LEE
    case "Alumina"
        spectral_band_edges = []; % For Alumina
end
N_band = max(length(spectral_band_edges)-1,1);

%% Define geometry
VS_TES = true(size_VS(1),size_VS(2),size_VS(3)); % Just a rectangle filling the domain
%Tz = linspace(T_TES_min,T_TES_i);
%for i =1:size_VS(3)
%    VS_T_init(:,:,i) = Tz(i);
%end
VS_T_init(VS_TES) = T_TES_i; % Initial temperature field in TES, everything outside the TES can just be 0 for now
VS_T_fixed = false(size_VS(1),size_VS(2),size_VS(3));
VS_T_fixed(~VS_TES) = true;

TES_system = TES_System(VS_TES,TES_material_name,vx_scale,VS_T_init);
TES_system.addFixedTemperature(VS_T_fixed);
TES_system.addDefaultBC();

%% Define our heater(s) in the TES
N_heater = N_heater_side^2;
loc_heaters = zeros(N_heater,2);
z_bot = size_VS(3)/2*(1-L_frac_heater);
z_top = size_VS(3)/2*(1+L_frac_heater);
r_vx_heater = r_heater/vx_scale(1); % voxel units
r_vx_void_heater = r_void_heater/vx_scale(1);
if N_heater_side > 1
    heater_x_spacing = size_VS(1)/(N_heater_side-1);
    heater_y_spacing = size_VS(2)/(N_heater_side-1);
else
    heater_x_spacing = size_VS(1)/2;
    heater_y_spacing = size_VS(2)/2;
end
for n1 = 1:N_heater_side
    if N_heater_side > 1
        x_loc = round(2*(n1-1)*heater_x_spacing)/2;
    else
        x_loc = size_VS(1)/2;
    end
    for n2 = 1:N_heater_side
        if N_heater_side > 1
            y_loc = round(2*(n2-1)*heater_y_spacing)/2;
        else
            y_loc = size_VS(2)/2;
        end
        TES_system.addHeater(x_loc,y_loc,z_bot,z_top,r_vx_heater,r_vx_void_heater,emissivity_heater,T_max_heater,0);
        nn = (n1-1)*N_heater_side+n2;
        loc_heaters(nn,:) = [x_loc,y_loc];
    end
end
heater_total_power = 0;
for nn = 1:N_heater
    heater_total_power = heater_total_power + TES_system.heaters{nn}.max_power;
end

%% Define our pipe(s) in the TES
air_table = readtable('air_properties.txt','Delimiter',' '); % Table of air properties vs T at p=1atm
N_pipe = N_pipe_side^2;
loc_pipes = zeros(N_pipe,2);
pipe_x_spacing = size_VS(1)/N_pipe_side;
pipe_y_spacing = size_VS(2)/N_pipe_side;

N_delete = 0;
for n1 =1:N_pipe_side
    for n2 = 1:N_pipe_side
        nn = (n1-1)*N_pipe_side+n2 - N_delete;
        loc_pipes(nn,:) = round(2*[(n1-0.5)*pipe_x_spacing,(n2-0.5)*pipe_y_spacing])/2;
        if ~ismember(loc_heaters,loc_pipes(nn,:),'rows')
            switch pipe_shape
                case "circle"
                    TES_system.addPipe(loc_pipes(nn,1),loc_pipes(nn,2),pipe_shape,r_vx_pipe,pipe_roughness,air_table);
                case "square"
                    TES_system.addPipe(loc_pipes(nn,1),loc_pipes(nn,2),pipe_shape,xy_width_vx_pipe,pipe_roughness,air_table);
            end
        else
            N_delete = N_delete + 1;
            loc_pipes(end,:) = [];
        end
    end
end


%% Define radiation domain
specular_BCs = true(2,3);
TES_system.createVoxelSpaces(spectral_band_edges,ns,specular_BCs);

%% Simulation

heat_capacity = TES_system.getHeatCapacity() % J/K
Delta_T_max = T_TES_max-T_TES_min;
energy_capacity = heat_capacity*Delta_T_max/1e3 % kJ
energy_density = energy_capacity/(xy_length_TES*xy_length_TES*z_height_TES)/1e3/3600 % MWh/m3
charge_rate_to_energy_capacity = heater_total_power/1e3/(energy_capacity/3600) % 1/h
nominal_discharge_power = nominal_discharge_power_fraction*energy_capacity/3600*1000 % W
cur_time = 0;

total_N_time_steps = sum(N_time_steps);
time = zeros(total_N_time_steps,1);
discharge_power = zeros(total_N_time_steps,1);
p_loss = zeros(total_N_time_steps,1);
m_dot = zeros(total_N_time_steps,1);
T_outlet = zeros(total_N_time_steps,1);
max_dT = zeros(total_N_time_steps,1);
count = 0;
tic_sim = tic;
power_drop_bool = false;
for nn = 1:length(N_time_steps)
    for i = 1:N_time_steps(nn)
        count = count+1;
        cur_time = cur_time + time_step_size(nn);
        if nn <= length(nominal_discharge_power)
            cur_discharge_power = nominal_discharge_power(nn);
        else
            cur_discharge_power = nominal_discharge_power(end);
        end
        if nn <= length(heater_power_fraction)
            cur_heater_power = heater_power_fraction(nn);
        else
            cur_heater_power = heater_power_fraction(end);
        end
        if nn <= length(N_rays)
            cur_N_rays = N_rays(nn);
        else
            cur_N_rays = N_rays(end);
        end


        tic_iter = tic;
        results = TES_system.chargeAndDischarge(time_step_size(nn),1,cur_N_rays,T_air_i,T_process,cur_discharge_power,cur_heater_power);
        disp(results)
        discharge_power(count) = -results.TES_power_rate;
        if discharge_power(count) < cur_discharge_power*0.99 && power_drop_bool == false
            cur_energy = heat_capacity*mean(TES_system.VS_T(TES_system.VS_TES_material)-T_TES_min,'all')/1e3;
            discharge_fraction_at_constant_power = (energy_capacity-cur_energy)/energy_capacity;
            power_drop_bool = true;
        end
        T_outlet(count) = results.T_fluid_out;
        p_loss(count) = results.pressure_loss;
        m_dot(count) = results.mass_flow_rate_per_pipe;
        time(count) = cur_time;
        time_iter = toc(tic_iter);
        if cur_time<300
            time_unit_fac = 1;
            time_unit = "seconds";

        else
            time_unit_fac = 60;
            time_unit = "minutes";
        end
        fprintf("current simulation time = %0.1f %s. Last step computer time = %0.1f seconds\n",cur_time/time_unit_fac,time_unit,time_iter)
        fprintf("steps remaining: %d\n",sum(N_time_steps)-count)
        if display_bool
            TES_system.visualizeTES(z_sample);
        end
        VS_T_cur = TES_system.VS_T;
        VS_T_cur(~TES_system.VS_TES_material) = nan;
        max_T_z = squeeze(max(VS_T_cur,[],[1,2],"omitnan"));
        min_T_z = squeeze(min(VS_T_cur,[],[1,2],"omitnan"));
        max_dT(count) = max(max_T_z-min_T_z);
    end
end
time_sim = toc(tic_sim)
final_results.simulation_time_mins = round((time_sim/60)*10)/10;
final_results.max_dT = round(max(max_dT)*10)/10;
final_results.discharge_fraction_at_constant_power = discharge_fraction_at_constant_power;
rho_in = 101325/287/T_air_i;
fan_power = p_loss.*m_dot*N_pipe/rho_in/0.8; % Assume efficiency of 80%

T_less_than_nom = find(T_outlet<T_process-273.15);

f = figure
tiledlayout(2,2)
nexttile(1);

plot(time/60,T_outlet);
hold on
yline(T_process-273.15,'--');
if ~isempty(T_less_than_nom)
    t_less_than_nom = time(T_less_than_nom(1));
    plot([t_less_than_nom,t_less_than_nom]/60,[0,T_process-273.15],'--','Color','k')
end
text((max(time)-min(time))/60*0.05,T_process-273.15+(max(T_outlet-273.15)*0.1),"T_{process}");
xlabel("t (min)")
ylabel("Outlet T (°C)")
ylim([0,inf])
nexttile(2);

plot(time/60,fan_power/1e3*N_module);
xlabel("t (min)")
ylabel("Fan power (kW)")
nexttile(3);

plot(time/60,discharge_power/1e6*N_module);
ylim([0,inf])
xlabel("t (min)")
ylabel("Discharge power (MW)")
nexttile(4)

plot(time/60,max_dT);
xlabel("t (min)")
ylabel("\Delta T (°C)")

count = 1;
figpath_base = fullfile("results","plots",sprintf("%s_plot",save_filename));
figpath = sprintf("%s_%02d.fig",figpath_base,count);
%while(isfile(figpath))
%    count = count+1;
%    figpath = sprintf("%s_%02d.fig",figpath,count);
%end
saveas(f,figpath,'fig')
save(fullfile("results","TES_objects",sprintf("%s_TES_obj_%d",save_filename,count)),"TES_system");

