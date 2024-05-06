%Transforms a direction vector from an hemispherical coordinate system
%built in a surface (with n normal to surface) to the right-hand GLOBAL
%coordinate system. 
% Written by Sebastian Sas Brunser with the help and input from Lukas
% Schlagenhauf 
% 20.05.2021


function [U_global]=transformCoord(U_local,n)

    %if zenith of hemisphere match the GLOBAL Z (^k), then no need to do anything:
    if (abs(n(1))<1e-5 && abs(n(2))<1e-5)
        if n(3)<0
            %R_y(180°)
            R=[-1 0 0;0 1 0;0 0 -1];
            U_global=(R*U_local').';
        else
            U_global=U_local;
        end
        return
        
    else

        nn=norm(n);
        nxy=sqrt(sum(n(1)^2+n(2)^2));

        cos_theta=n(3)/nn;
        sin_theta=nxy/nn;

        cos_phi=n(1)/nxy;
        sin_phi=n(2)/nxy;


        %Rotation Matrix, defined from elemental rotation matrices:
        %R=inv(R_y(-theta)*R_z(-phi))
        
        R=[cos_theta*cos_phi    -sin_phi       sin_theta*cos_phi;...
          cos_theta*sin_phi     cos_phi         sin_theta*sin_phi;...
          -sin_theta            0           cos_theta];


        U_global=(R*U_local').';
        
    end
end