clear
close all
clc

%% USER INPUT PARAMETERS
% Geometry
N_grid_D = 24; % Number of gridcells across the diameter
width_TES = 0.1; % [m]: x-direction
length_TES = 0.1; % [m]: y-direction
height_TES = 0.1; % [m]: z-direction

% Discharge pipe parameters
% Number of pipes in the each direction in the domain (pipe axis always in z-direction)
% pipes are automatically generated, evenly spaced based on N_pipe_side
% If the automatically generated locations results in a pipe and heater at the same location, the pipe is not generated
N_pipe_side = 1;
r_pipe = 0.1*width_TES;
pipe_roughness = 0.0015/1000; % m: 0.0015 for glass, source: https://www.engineersedge.com/fluid_flow/pipe-roughness.htm
T_air_i = 20 + 273.15; % Kelvin

% TES Parameters - borosilicate glass SHOTT-BK7 - This is now handled by TES_material class!
TES_material_name = "SHOTT_BK7_LEE"; % SHOTT_BK7_LEE, JGS2
%TES_material_name = "JGS2"

N_rays = 2e6;
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
T_TES_i = 100 + 273.15; % Kelvin: Initial temperature of the TES
N_time_steps = [50];
time_step_size = [50]; % seconds
heater_power_fraction = [0];
nominal_discharge_power_fraction = [1/5]; % 1/h: Fraction of storage capacity
T_process = 500 + 273.15; % Kelvin
z_sample = 1;

%%%%% END USER DEFINED PARAMETERS %%%%%
%% Check inputs
if N_heater_side == N_pipe_side
    error("N_heater_" + ...
        "side cannot equal N_pipe_side")
end

%% Derived parameters
vx_scale = width_TES/N_grid_D; % side length of each grid cell
N_grid_z = round(height_TES/vx_scale);
height_TES = vx_scale*N_grid_z; % ,m, Actual L;
r_vx_pipe = r_pipe/vx_scale;
D_pipe = r_pipe*2;

void_fraction = (pi*r_pipe^2)/(width_TES*length_TES);
volumetric_fraction = 1-void_fraction

VS_T_init = zeros(N_grid_D,N_grid_D,N_grid_z);
size_VS = size(VS_T_init);


switch TES_material_name
    case "JGS2"
        spectral_band_edges = [0.22,0.27,2,2.5,3,3.5,4.5,1000]; % good band edges for JGS2
    case "SHOTT_BK7_LEE"
        spectral_band_edges = [0.38,1.6,2.03,2.68,1000]; % Good band edges for SHOTT BK7 LEE
    case "Alumina"
        spectral_band_edges = []; % For Alumina
end
N_band = max(length(spectral_band_edges)-1,1);

%% Define geometry
VS_TES = true(size_VS); % Just a rectangle filling the domain
VS_T_init(VS_TES) = T_TES_i; % Initial temperature field in TES, everything outside the TES can just be 0 for now
VS_T_fixed = false(size_VS);
VS_T_fixed(~VS_TES) = true;

TES_system = TES_System(VS_TES,TES_material_name,vx_scale,VS_T_init);
TES_system.addFixedTemperature(VS_T_fixed);
TES_system.addDefaultBC();

%% Define our heater(s) in the TES
N_heater = N_heater_side^2;
loc_heaters = zeros(N_heater,2);
z_bot = size_VS(3)/2*(1-L_frac_heater);
z_top = size_VS(3)/2*(1+L_frac_heater);
r_vx_heater = r_heater/vx_scale; % voxel units
r_vx_void_heater = r_void_heater/vx_scale;
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
            TES_system.addPipe(loc_pipes(nn,1),loc_pipes(nn,2),r_vx_pipe,pipe_roughness,air_table);
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
Delta_T_max = TES_system.TES_material.annealing_temperature-T_TES_min
energy_capacity = heat_capacity*Delta_T_max/1e3 % kJ
energy_density = energy_capacity/(width_TES*length_TES*height_TES)/1e3/3600 % MWh/m3
charge_rate_to_energy_capacity = heater_total_power/1e3/(energy_capacity/3600) % 1/h
nominal_discharge_power = nominal_discharge_power_fraction*energy_capacity/3600*1000 % W
cur_time = 0;

total_N_time_steps = sum(N_time_steps);
time = zeros(total_N_time_steps,1);
discharge_power = zeros(total_N_time_steps,1);
p_loss = zeros(total_N_time_steps,1);
m_dot = zeros(total_N_time_steps,1);
outlet_T = zeros(total_N_time_steps,1);
max_dT = zeros(total_N_time_steps,1);
count = 0;
for nn = 1:length(N_time_steps)
    for i = 1:N_time_steps(nn)
        count = count+1;
        cur_time = cur_time + time_step_size(nn);
        
        results = TES_system.chargeAndDischarge(time_step_size(nn),1,N_rays,T_air_i,T_process,nominal_discharge_power(nn),heater_power_fraction(nn));
        disp(results)
        discharge_power(count) = results.TES_power_rate;
        outlet_T(count) = results.T_fluid_out;
        p_loss(count) = results.pressure_loss;
        m_dot(count) = results.mass_flow_rate_per_pipe;
        time(count) = cur_time;
        fprintf("Current time = %0.0f seconds\n",cur_time)
        TES_system.visualizeTES(z_sample);
        max_T_z = squeeze(max(TES_system.VS_T,[],[1,2]));
        min_T_z = squeeze(min(TES_system.VS_T,[],[1,2]));
        max_dT(count) = max(max_T_z-min_T_z);
    end
end
rho_in = 101325/287/T_air_i;
fan_power = p_loss.*m_dot*N_pipe/rho_in/0.8; % Assume efficiency of 80%
figure
tiledlayout(2,2)
nexttile(1);
plot(time/60,outlet_T);
nexttile(2);
plot(time/60,-discharge_power);
ylim([0,inf])
nexttile(3);
plot(time/60,fan_power);
nexttile(4)
plot(time/60,max_dT);