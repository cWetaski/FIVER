function [stress_comps,vm_stress] = getThermalStressCylinder(r,a,b,Delta_T,alpha,nu,E)
%GETTHERMALSTRESSCYLINDER Summary of this function goes here
%   Delta T = T_a - T_b
sigma_r = alpha*E*Delta_T/(2*(1-nu)*log(b/a))*(-log(b/r)+(a^2*(b^2-r^2))/(r^2*(b^2-a^2))*log(b/a));
sigma_theta = alpha*E*Delta_T/(2*(1-nu)*log(b/a))*(1-log(b/r)-(a^2*(b^2+r^2))/(r^2*(b^2-a^2))*log(b/a));
sigma_z = sigma_r + sigma_theta;
stress_comps = [sigma_r(:),sigma_theta(:),sigma_z(:)];
vm_stress = sqrt(1/2*((sigma_r-sigma_theta).^2+(sigma_r-sigma_z).^2+(sigma_theta-sigma_z).^2));
end

