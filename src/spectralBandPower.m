%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07 (Charles Wetaski)

function E_b_band = spectralBandPower(lambda1,lambda2,T,nn)
    % SPECTRALBANDPOWER
    % Integrates spectral blackbody emissive power between two specified
    % wavelengths in um for a given T and refractive index nn
    % Inputs:
    %   lambda1 (scalar double) [um]:   Lower wavelength of spectral band
    %   lambda2 (scalar double) [um]:   Upper wavelength of spectral band (lambda1 <= lambda2)
    %   T (double) [K]:                 Temperature in kelvin. Can be scalar or an array of any dimension, output will have
    %                                   same shape as T.
    %   nn (double) [-]:                Refractive index of medium (or medium the surface is emitting into) in the band
    %                                   Can be scalar or array in same dimension of T
    % Outputs:
    %   E_b_band: [W]:                  Blackbody emissive power in the specified wavelength band
    
    %% Constants
    h = 6.62607004*10^(-34);     % Planck's constant [J s]
    c = 299792458;              % speed of light in a vacuum [m/s]
    Kb = 1.38064852*10^(-23);    % Boltzmann's constant [J/K]
    lambda1_m = lambda1*10^(-6); % wavelength in meters
    lambda2_m = lambda2*10^(-6);
    
    
    %% Planck's Law of EM Radiation [W/(m^2 um)]
    f = @(lambda_m) (2*pi*h*c^2)./(nn.^2*lambda_m.^5.*(exp((h*c)./(nn.*lambda_m.*Kb.*T))-1)); % [W/(m^2-m)]
    E_b_band = integral(f,lambda1_m,lambda2_m,'ArrayValued',true);
end


