% ©2024 ETH Zurich, Charles Wetaski, Sebastian Sas Brunser, Emiliano Casati
%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-09-06 (Charles Wetaski)
% 

function Q_emit = totalEmissivePowerSingle(spectral_band_edges,voxel_spaces,T,ind)
    %TOTALEMISSIVEPOWERSINGLE Get total emissive power at a single voxel in a spectral voxel space at a given temperature
    %   inputs:
    %       spectral_band_edges (1D double) [um]:   Edges of spectral bands (in increasing wavelength)
    %       voxel_spaces (1D cell):                 Cell array of voxel_space objects in each spectral band
    %       T (scalar double)    [K]:               Temperature to evaluate the emissive power at
    %       ind (scalar double (int)):              Linear index of the voxel of interest.
    %   outputs:
    %       Q_emit (scalar) [W]: Total emissive power of the voxel  
Q_emit = 0;
if T == 0
    return;
end
N_bands = length(spectral_band_edges)-1;
vx_scale = voxel_spaces{1}.voxel_scale;
for n = 1:N_bands
    nn = voxel_spaces{n}.refractive_indexes(ind);
    lb = spectral_band_edges(n-1);
    ub = spectral_band_edges(n);
    PM_kappa = voxel_spaces{n}.PM_absorption_coeffs(ind);
    if PM_kappa == 0
        surf_area = voxel_spaces{n}.surface_areas(ind);
        if surf_area == 0
            continue
        else
            opaq_eps = voxel_spaces{n}.opaque_emissivities(ind);
            if opaq_eps == 0
                continue;
            else
            Q_emit = Q_emit + vx_scale.^2*surf_area.*opaq_eps.*spectralBandPower(lb,ub,T,nn);
            
            end
        end
    else
        Q_emit = Q_emit + vx_scale.^2*4*PM_kappa.*spectralBandPower(lb,ub,T,nn);
    end
end

