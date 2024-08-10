function [M_pipe, RHS_BC_pipe] = PipeBC_1D(mesh,ii_pipe,k,h,T_b)
% PIPEBC1D Define a convective boundary condition at index ii which gives the heat transfer equivalent to a cylindrical
% pipe (its axis orthogonal to the x-axis) which has heat transfer coefficient h through a solid with heat conductivity k
% It is assumed that flow is symmetric along the perimiter of the cylinder and heat conductivity is the same in both
% boundary cells. Note that this assumes the uniform cartesian grid!
% This code is a direct rewriting of the boundaryCondition1D code from FVTool, but we are essentially have a BC inside
% the domain!
% Credit: Eftekhari, A.A. et al. (2015). FVTool: a finite volume toolbox for Matlab. Zenodo. http://doi.org/10.5281/zenodo.32745
% Parameters:
%   mesh: Domain mesh variable
%   T: T_cell: (1 x Nx) Cell variable of temperature distribution in system before creating this internal BC
%   ii: i-index of pipe location
%   jj: j-index of pipe location
%   k [W/(K-m)]: (scalar) Heat conductivity of the wall material
%   h [W/(K-m^2))]: (scalar) Heat transfer coefficient fluid at each z-location along the pipe axis
%   T_b [K]: (scalar) Fluid temperature of fluid at pipe location.

Nxyz = mesh.dims;
Nx = Nxyz(1);
G= 1:(Nx+2);
dx = mesh.cellsize.x(ii_pipe+1);

nb = 3; % 3 nodes (2 boundary nodes and the center node)

% define the vectors to be used for the creation of the sparse matrix
iii = zeros(nb,1);
jjj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
RHS_BC_pipe = zeros(Nx+2,1);

% These values obtained from reformulating a boundary condition of the third kind into a*dT/dx+b*T=c
% NOTE!!! when the surface normal is negative (i.e., in 1D, when the surface is located on the right side of the
% boundary), we use -a instead of a!!!!
% a = k;
% b = h*pi/4; % Since the CFD solver is using rectangular cartesian cell boundaries, we modifiy the heat transfer coefficient by multiplying it by the ratio of perimiters
%              % since the heat transfer coefficient is determined for a cylindrical pipe with diameter equal to the cellsize 
% c = h*pi/4*T_b; %
% 
% K1 = a/dx+b/2; % Define these constants to simplify algebra (which I may have to type out formally somewhere!!) Note, (b = K1+K2)
% K2 = -a/dx+b/2;

h_star = h*pi/4;

K1 = -k/dx;
K2 = k/dx - h_star*(1-log(2)/(2*k/h/dx+log(2))); % Note regular h in denominator, not h_star!
K3 = -h_star*T_b*(1-log(2)/(2*k/h/dx+log(2)));


ii = ii_pipe+1; % pipe location in cell array (i.e., including the ghost cells, we have to shift+1)
% Define matrix of coefficients and RHS for ghost cell for solving syste

% At pipe ('ghost cell')
q = 1;
iii(q) = G(ii);  jjj(q) = G(ii); s(q) = K1; % at pipe
RHS_BC_pipe(G(ii)) = K3;
% left boundary
q=q+1;
iii(q) = G(ii);  jjj(q) = G(ii-1);  s(q) = K2/2; % left of pipe
% right boundary
q=q+1;
iii(q) = G(ii); jjj(q) = G(ii+1); s(q) = K2/2; % right of pipe


% Build the sparse matrix of the boundary conditions
M_pipe = sparse(iii(1:q), jjj(1:q), s(1:q), Nx+2, Nx+2);
end

