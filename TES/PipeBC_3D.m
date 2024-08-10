function [M_pipe, RHS_BC_pipe] = PipeBC_3D(mesh,ii_pipe,jj_pipe,k,h,T_b)
% PIPEBC1D Define a convective boundary condition at index ii which gives the heat transfer equivalent to a cylindrical
% pipe (its axis orthogonal to the x-axis) which has heat transfer coefficient h through a solid with heat conductivity k
% It is assumed that flow is symmetric along the perimiter of the cylinder and heat conductivity is the same in both
% boundary cells. Note that this assumes the uniform cartesian grid in the x and y axis!
% This code is a direct rewriting of the boundaryCondition1D code from FVTool, but we are essentially have a BC inside
% the domain!
% Credit: Eftekhari, A.A. et al. (2015). FVTool: a finite volume toolbox for Matlab. Zenodo. http://doi.org/10.5281/zenodo.32745
% Parameters:
%   mesh: Domain mesh variable
%   ii: i-index of pipe location
%   jj: j-index of pipe location
%   k [W/(K-m)]: (1 x Nz) Heat conductivity of the wall material
%   h [W/(K-m^2))]: (scalar) Heat transfer coefficient fluid at each z-location along the pipe axis
%   T_b [K]: (scalar) Fluid temperature of fluid at pipe location.

Nxyz = mesh.dims;
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2,Nz+2);
dx = mesh.cellsize.x(ii_pipe+1);
dy = mesh.cellsize.y(jj_pipe+1);
if dx ~= dy
    disp("Warning, nonuninform mesh in xy-plane (PipeBC_3D)")
end
dxy = dx;

% number of nodes
nb = 5*Nz; % (1 center node and 4 boundary nodes for each z-point)

% define the vectors to be used for the creation of the sparse matrix
iii = zeros(nb,1);
jjj = zeros(nb,1);
s = zeros(nb,1);

RHS_BC_pipe = zeros((Nx+2)*(Ny+2)*(Nz+2),1);

% % These values obtained from reformulating a boundary condition of the third kind into a*dT/dx+b*T=c
% a = -k;
% b = -h*pi/4; % Since the CFD solver is using rectangular cartesian cell boundaries, we modifiy the heat transfer coefficient by multiplying it by the ratio of perimiters
%              % since the heat transfer coefficient is determined for a cylindrical pipe with diameter equal to the cellsize 
% c = -h*pi/4.*T_b; %
% 
% K1 = a/dxy+b/2; % Define these constants to simplify algebra (which I may have to type out formally somewhere!!) Note, (b = K1+K2)
% K2 = -a/dxy+b/2;

% NEW VERSION!
h_star = h*pi/4;

K1 = -k/dxy;
K2 = k/dxy - h_star.*(1-log(2)./(2*k./h/dxy+log(2))); % Note regular h in denominator, not h_star!
K3 = -h_star.*T_b.*(1-log(2)./(2*k./h/dxy+log(2)));

ii = ii_pipe+1; % pipe location in cell array (i.e., including the ghost cells, we have to shift+1)
jj = jj_pipe+1;

% Define matrix of coefficients and RHS for ghost cell for solving system

% Obtain the columns we want in advance so we don't need to index the 3D matrix more than necessary
iii_pipe = G(ii,jj,2:end-1); iii_left = G(ii-1,jj,2:end-1); iii_right = G(ii+1,jj,2:end-1); iii_bot = G(ii,jj-1,2:end-1); iii_top = G(ii,jj+1,2:end-1); iii_back = G(ii,jj,2:end-2); iii_front = G(ii,jj,3:end-1); 

% at pipe
q = 1:Nz;
iii(q) = iii_pipe; jjj(q) = iii_pipe; s(q) = K1; % at pipe
RHS_BC_pipe(iii_pipe) = K3;
% left boundary
q=q(end)+(1:Nz);
iii(q) = iii_pipe; jjj(q) = iii_left;  s(q) = K2/4; % left of pipe. 
% right boundary
q=q(end)+(1:Nz);
iii(q) = iii_pipe; jjj(q) = iii_right; s(q) = K2/4; % right of pipe
% bottom boundary
q = q(end)+(1:Nz);
iii(q) = iii_pipe; jjj(q) = iii_bot;  s(q) = K2/4; % bottom of pipe.
% top boundary
q=q(end)+(1:Nz);
iii(q) = iii_pipe; jjj(q) = iii_top; s(q) = K2/4; % top of pipe

% Build the sparse matrix of the boundary conditions
M_pipe = sparse(iii, jjj, s, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));

end

