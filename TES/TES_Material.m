classdef TES_Material
    %TES_Material This class stores a bunch of properties for candidate TES materials
    %   The materials don't all have properties in the same form. e.g., index of refraction is given for
    %   SHOTT BK7 by the Sellmeier formula, but by interpolation from a data table for quartz.
    
    properties
        material_name
        thermal_conductivity        % W/m-K
        density                     % kg/m3
        specific_heat               % J/(kg-K)
        annealing_temperature       % K
        softening_temperature       % K
        emissivity_coefs            % Coefficients for emissivity formula
        spectral_absorptivity_data  % This should be a table
        refraction_mode;            % Table/Sellmeier
        refractive_index_data       % This should be a table, alternatively use Sellmeier equation (but not both!)
        dispersion_formula_B        % [B1,B2,B3]: Sellmeier equation
        dispersion_formula_C        % [C1,C2,C3] um^2: Sellmeir equation, 
    end
    
    methods
        function obj = TES_Material(material_name)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            materials = TES_Material.getAvailableMaterials;
            if ~ismember(material_name,materials)
                error("Specified material name not available,\nAvailable materials: %s",strjoin(materials,", "))
            end

            obj.material_name = material_name;
            
            switch material_name
                case 'SHOTT_BK7_LEE' % https://www.schott.com/shop/advanced-optics/en/Optical-Glass/SCHOTT-N-BK7/c/glass-SCHOTT%20N-BK7%C2%AE
                    obj.thermal_conductivity = 1.114; % W/(K-m); Shott BK7
                    obj.density = 2510;
                    obj.specific_heat = 858;
                    obj.annealing_temperature = 557+273.15;
                    obj.softening_temperature = 719+273.15;
                    obj.dispersion_formula_B = [1.03961212,0.231792344,1.01046945];
                    obj.dispersion_formula_C = [0.00600069867,0.0200179144,103.560653];
                    obj.refraction_mode = "Sellmeier";
                    obj.spectral_absorptivity_data = readtable("BK7_data_Lee.txt"); % Note: BK7_data_Lee comes from Lee and Viskanta, 1998, Fig. 2, and gives values in 1/cm: BK7 data from SHOTT is internal transmittance and only goes to 2.5 um
                case 'SHOTT_BK7'% https://www.schott.com/shop/advanced-optics/en/Optical-Glass/SCHOTT-N-BK7/c/glass-SCHOTT%20N-BK7%C2%AE
                    obj.thermal_conductivity = 1.114; % W/(K-m); Shott BK7
                    obj.density = 2510;
                    obj.specific_heat = 858;
                    obj.annealing_temperature = 557+273.15;
                    obj.softening_temperature = 719+273.15;
                    obj.dispersion_formula_B = [1.03961212,0.231792344,1.01046945];
                    obj.dispersion_formula_C = [0.00600069867,0.0200179144,103.560653];
                    obj.refraction_mode = "Sellmeier";
                    obj.spectral_absorptivity_data = readtable("BK7_data.txt"); % Note: BK7_data_Lee comes from Lee and Viskanta, 1998, Fig. 2, and gives values in 1/cm: BK7 data from SHOTT is internal transmittance and only goes to 2.5 um
                case 'JGS2'
                    obj.thermal_conductivity = 1.4; % W/(K-m);
                    obj.density = 2200;
                    obj.specific_heat = 670;
                    obj.annealing_temperature = 1215+273.15;
                    obj.softening_temperature = 1683+273.15;
                    obj.refractive_index_data = readtable("JGS2_refraction.txt");
                    obj.refraction_mode = "Table";
                    obj.spectral_absorptivity_data = readtable("JGS2_transmittance.txt"); % Note: BK7_data_Lee comes from Lee and Viskanta, 1998, Fig. 2, and gives values in 1/cm: BK7 data from SHOTT is internal transmittance and only goes to 2.5 um       
                case 'Alumina'
                    obj.thermal_conductivity = 5.7; % 37 @ 25°C, 10 @ 500C, 5.7 at 1200C
                    obj.density = 3987;
                    obj.specific_heat = 840;
                    obj.annealing_temperature=2054+273.15;
                    obj.softening_temperature=2054+273.15;
                    obj.emissivity_coefs = [0.98, -53, 10.2]; % Source: Jones (2018), A COMPILATION OF DATA ON THE RADIANT EMISSIVITY OF SOME MATERIALS AT HIGH TEMPERATURES
                end
        end
        function mean_abs_coeff = getAbsorptivity(obj,wavelength1,wavelength2)
            if strcmp(obj.material_name,'Alumina')
                mean_abs_coeff = -1;
                return
            end
            data_table = obj.spectral_absorptivity_data;
            if wavelength1 < data_table.wavelength(1) || wavelength1 > data_table.wavelength(end)
                mean_abs_coeff = -1; % Cannot interpolate
            else
                if strcmp(obj.material_name,'SHOTT_BK7_LEE') % Data in 1/cm
                    fun = @(x) interp1(data_table.wavelength,data_table.abs_coeff,x)*100; % 1/m
                elseif strcmp(obj.material_name,'SHOTT_BK7') % Data as extinction coefficient
                    fun = @(x) 4*pi*interp1(data_table.wavelength,data_table.k,x)/(x/1e6); % 1/m
                elseif strcmp(obj.material_name,'JGS2') % Data as transmittivity (%) (ASSUMING 10mm SAMPLE FOR NOW!!!)
                    fun = @(x) -log(interp1(data_table.wavelength,data_table.transmissivity,x)/100)/(0.01); % 1/m
                else % If using a different data source: must add to this line!
                    error("NO DATA SOURCE")
                end
                if nargin == 2 % Just return the value 
                    mean_abs_coeff = fun(wavelength1);
                else % Take mean in wavelength band
                    if wavelength2 < data_table.wavelength(1) || wavelength2 > data_table.wavelength(end) || wavelength2 <= wavelength1
                        mean_abs_coeff = -1; % Cannot interpolate
                    else
                        mean_abs_coeff = integral(fun,wavelength1,wavelength2)/(wavelength2-wavelength1);
                    end
                end
            end
        end
        
        function n = getRefractiveIndex(obj,wavelength1,wavelength2)
            if strcmp(obj.material_name,"Alumina")
                n = 1;
                return
            end
            if strcmp(obj.refraction_mode,"Sellmeier")
                B = obj.dispersion_formula_B;
                C = obj.dispersion_formula_C;
                
                fun = @(x) real(sqrt(1 + B(1)*x.^2./(x.^2-C(1)) ...
                                       + B(2)*x.^2./(x.^2-C(2)) ...
                                       + B(3)*x.^2./(x.^2-C(3))));
            elseif strcmp(obj.refraction_mode,"Table")
                data_table = obj.refractive_index_data;
                fun = @(x) interp1(data_table.wavelength,data_table.n,x);
            end
            if nargin == 2
                n = fun(wavelength1);
            else
                n = integral(fun,wavelength1,wavelength2)/(wavelength2-wavelength1);
            end
        end
        function epsilon = getEmissivity(obj,T)
            if strcmp(obj.material_name,"Alumina")
                epsilon = obj.emissivity_coefs(1)+obj.emissivity_coefs(2)*T*1e-5+obj.emissivity_coefs(3)*T.^2*1e-8; % Note, this doesn't account for any spectral behaviour!
            else
                epsilon = 1;
            end
        end

    end
    methods(Static)
        function materials = getAvailableMaterials()
            materials = ["SHOTT_BK7_LEE","SHOTT_BK7","JGS2","Alumina"];
        end

        function listAvailableMaterials() % This function must be manually updated!
            disp(TES_Material.getAvailableMaterials())
        end
    end
end