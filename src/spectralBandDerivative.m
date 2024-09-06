% ©2024 ETH Zurich, Charles Wetaski, Sebastian Sas Brunser, Emiliano Casati
%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-09-06 (Charles Wetaski)

function dEdT = spectralBandDerivative(lambda1,lambda2,T,dT)
    % SPECTRALBANDDERIVATIVE
    % Computes derivative (with respect to temperature) of spectral blackbody emissive power between two specified
    % wavelengths in um at a given T
    % Inputs:
    %   lambda1 (scalar double) [um]:   Lower wavelength of spectral band
    %   lambda2 (scalar double) [um]:   Upper wavelength of spectral band (lambda1 <= lambda2)
    %   T (double) [K]:                 Temperature in kelvin. Can be scalar or an array of any dimension, output will have
    %                                   same shape as T.
    % (optional)
    %   dT (double) [K]:                Temperature difference over which to evaluate the derivative, otherwise determined automatically.
    %
    % Outputs:
    %   dE_b_band: [W]:                  Blackbody emissive power in the specified wavelength band
    
    %% Constants
    if nargin == 3
        dT = T*sqrt(eps); % eps here is machine epsilon, apparently this optimal due to some floating point stuff!
        dT(T==0) = sqrt(eps);
    end
    f_l = spectralBandPower(lambda1,lambda2,T,nn);
    f_u = spectralBandPower(lambda1,lambda2,T+dT,nn);
    dEdT = (f_u-f_l)./(dT);
end


