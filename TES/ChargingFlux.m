classdef ChargingFlux < InternalFlux 
    %CHARGINGFLUX A subclass of Internal Flux which can easily be assigned and removed from spectral voxel space
    %   Importantly, the spectral behaviour of the flux is a function of the power and surface area, so the power can be adjusted
    %   and the spectral behaviour will automatically change.
    %   A chargingflux is a associated with a parent "Heater" object, which is used to assign it to a voxel space and remove it depending on
    %   whether the TES system is charging.
    
    properties
        heater; % parent heater object
    end
    
    methods
        function obj = ChargingFlux(power,ray_gen_function,heater)
            %CHARGINGFLUX Construct an instance of this class
            %   Detailed explanation goes here
            obj@InternalFlux(power,ray_gen_function);
            obj.heater = heater;
            
        end
    end
end

