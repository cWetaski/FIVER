function Q_emit = totalEmissivePowerSingle(lambdas,voxel_spaces,T,i)
%TOTALEMISSIVEPOWER Get total emissive power for voxel space with varying spectral properties
%   Detailed explanation goes here
Q_emit = 0;
if T == 0
    return;
end
N_bands = length(lambdas)+1;
vx_scale = voxel_spaces{1}.voxel_scale;
for n = 1:N_bands
    nn = voxel_spaces{n}.refractive_indexes(i);
    if n == 1
        lb = 0;
        ub = lambdas(n);
    elseif n == N_bands
        lb = lambdas(n-1);
        ub = 10^12/T/nn; % in case nan, just set it to 10^12
    else
        lb = lambdas(n-1);
        ub = lambdas(n);
    end
    PM_kappa = voxel_spaces{n}.PM_absorption_coeffs(i);
    if PM_kappa == 0
        surf_area = voxel_spaces{n}.surface_areas(i);
        if surf_area == 0
            continue
        else
            opaq_eps = voxel_spaces{n}.opaque_emissivities(i);
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

