%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07 (Charles Wetaski)

function dEdT = spectralBandDerivative(lambda1,lambda2,T,nn,dT)
    % SPECTRALBANDDERIVATIVE
    % Computes derivative (with respect to temperature) of spectral blackbody emissive power between two specified
    % wavelengths in um at a given T and refractive index nn
    % Inputs:
    %   lambda1 (scalar double) [um]:   Lower wavelength of spectral band
    %   lambda2 (scalar double) [um]:   Upper wavelength of spectral band (lambda1 <= lambda2)
    %   T (double) [K]:                 Temperature in kelvin. Can be scalar or an array of any dimension, output will have
    %                                   same shape as T.
    %   nn (double) [-]:                Refractive index of medium (or medium the surface is emitting into) in the band
    %                                   Can be scalar or array in same dimension of T
    % Outputs:
    %   dE_b_band: [W]:                  Blackbody emissive power in the specified wavelength band
    
    %% Constants
    if nargin == 4
        dT = T*sqrt(eps);
        dT(T==0) = sqrt(eps);
    end
    f_l = spectralBandPower(lambda1,lambda2,T,nn);
    f_u = spectralBandPower(lambda1,lambda2,T+dT,nn);
    dEdT = (f_u-f_l)./(dT);
end


