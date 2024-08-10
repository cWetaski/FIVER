function [M_pipe, RHS_BC_pipe] = PipeBC_2D(mesh,ii_pipe,jj_pipe,k,h,T_b)
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
Nx = Nxyz(1); Ny = Nxyz(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
dx = mesh.cellsize.x(ii_pipe+1);
dy = mesh.cellsize.y(jj_pipe+1);
if dx ~= dy
    disp("Warning, nonuninform mesh (PipeBC_2D)")
end
dxy = dx;

% number of boundary nodes (axact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 5; % 5 boundary nodes (1 center node and 4 boundary nodes)

% define the vectors to be used for the creation of the sparse matrix
iii = zeros(nb,1);
jjj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
RHS_BC_pipe = zeros((Nx+2)*(Ny+2),1);

% NEW VERSION!
h_star = h*pi/4;

K1 = -k/dxy;
K2 = k/dxy - h_star*(1-log(2)/(2*k./h/dxy+log(2))); % Note regular h in denominator, not h_star!
K3 = -h_star*T_b*(1-log(2)/(2*k./h/dxy+log(2)));

ii = ii_pipe+1; % pipe location in cell array (i.e., including the ghost cells, we have to shift+1)
jj = jj_pipe+1;
% Define matrix of coefficients and RHS for ghost cell for solving system
q = 1;
% at pipe
iii(q) = G(ii,jj);  jjj(q) = G(ii,jj); s(q) = K1; % at pipe
RHS_BC_pipe(G(ii,jj)) = K3;
% left boundary
q=q+1;
iii(q) = G(ii,jj);  jjj(q) = G(ii-1,jj);  s(q) = K2/4; % left of pipe. 
% right boundary
q=q+1;
iii(q) = G(ii,jj); jjj(q) = G(ii+1,jj); s(q) = K2/4; % right of pipe
% bottom boundary
q = q+1;
iii(q) = G(ii,jj);  jjj(q) = G(ii,jj-1);  s(q) = K2/4; % bottom of pipe.
% top boundary
q=q+1;
iii(q) = G(ii,jj); jjj(q) = G(ii,jj+1); s(q) = K2/4; % top of pipe


% Build the sparse matrix of the boundary conditions
M_pipe = sparse(iii(1:q), jjj(1:q), s(1:q), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));


end

