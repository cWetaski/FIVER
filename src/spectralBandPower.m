% ©2024 ETH Zurich, Charles Wetaski, Sebastian Sas Brunser, Emiliano Casati
%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-09-06 (Charles Wetaski)
%   
function E_b_band = spectralBandPower(lambda1,lambda2,T)
    % SPECTRALBANDPOWER
    % Integrates spectral blackbody emissive power between two specified
    % wavelengths in um for a given T for refractive index = 1 (vacuum)
    % We only consider the vacuum emissive power since spectral properties of real materials are typically reported
    % with respect to the vacuum wavelength, not the wavelength within the medium. Increases in the spectral band power
    % due to a nonunit refractive index (i.e., multiplication by n^2) must applied outside of this function.
    % Inputs:
    %   lambda1 (scalar double) [um]:   Lower wavelength of spectral band
    %   lambda2 (scalar double) [um]:   Upper wavelength of spectral band (lambda1 <= lambda2)
    %   T (double) [K]:                 Temperature in kelvin. Can be scalar or an array of any dimension, output will have
    %                                   same shape as T.
    % Outputs:
    %   E_b_band: [W]:                  Blackbody emissive power in the specified wavelength band
    
    %% Constants
    h = 6.62607004*10^(-34);     % Planck's constant [J s]
    c = 299792458;              % speed of light in a vacuum [m/s]
    Kb = 1.38064852*10^(-23);    % Boltzmann's constant [J/K]
    lambda1_m = lambda1*10^(-6); % wavelength in meters
    lambda2_m = lambda2*10^(-6);
    
    
    %% Planck's Law of EM Radiation [W/(m^2 um)]
    f = @(lambda_m) (2*pi*h*c^2)./(lambda_m.^5.*(exp((h*c)./(lambda_m.*Kb.*T))-1)); % [W/(m^2-m)]
    E_b_band = integral(f,lambda1_m,lambda2_m,'ArrayValued',true);
end


