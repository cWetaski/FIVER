function [M_pipe, RHS_BC_pipe] = PipeBC_2D_Multi(mesh,x_center,y_center,ri,k,h,T_b)
% PIPEBC1D Define a convective boundary condition at index ii which gives the heat transfer equivalent to a cylindrical
% pipe (its axis orthogonal to the x-axis) which has heat transfer coefficient h through a solid with heat conductivity k
% It is assumed that flow is symmetric along the perimiter of the cylinder and heat conductivity is the same in both
% boundary cells. Note that this assumes the uniform cartesian grid!
% This code is a direct rewriting of the boundaryCondition1D code from FVTool, but we are essentially have a BC inside
% the domain!
% Credit: Eftekhari, A.A. et al. (2015). FVTool: a finite volume toolbox for Matlab. Zenodo. http://doi.org/10.5281/zenodo.32745
% Parameters:
%   mesh: Domain mesh variable
%   ii: i-index of pipe location
%   jj: j-index of pipe location
%   k [W/(K-m)]: (scalar) Heat conductivity of the wall material
%   h [W/(K-m^2))]: (scalar) Heat transfer coefficient fluid at each z-location along the pipe axis
%   T_b [K]: (scalar) Fluid temperature of fluid at pipe location.

[VS_pipe,VS_boundaries,inds_pipe_boundary,inds_pipe_interior,boundary_perimiter_span] = getMultiPipe(mesh.dims,x_center,y_center,ri);
Nxyz = mesh.dims;
Nx = Nxyz(1); Ny = Nxyz(2);
ii_center = floor(x_center)+1;
jj_center = floor(y_center)+1;


G = reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
dx = mesh.cellsize.x(ii_center+1);
dy = mesh.cellsize.y(jj_center+1);

if dx ~= dy
    disp("Warning, nonuninform mesh (PipeBC_2D)")
end
dxy = dx;

% define the vectors to be used for the creation of the sparse matrix
N_boundary = size(inds_pipe_boundary,1);
N_interior = size(inds_pipe_interior,1);
nb = (N_boundary*5 + N_interior);
iii = zeros(nb,1);
jjj = zeros(nb,1);
s = zeros(nb,1);

% define the RHS column vector
RHS_BC_pipe = zeros((Nx+2)*(Ny+2),1);

% Define matrix of coefficients and RHS for ghost cell for solving system
q = 0;
% interior cells -> have no information set temperature to 0;
for nn = 1:N_interior
    
    ii = inds_pipe_interior(nn,1)+1;
    jj = inds_pipe_interior(nn,2)+1;
    q=q+1;
    iii(q) = G(ii,jj);  jjj(q) = G(ii,jj); s(q) = 1; % at pipe 
    RHS_BC_pipe(G(ii,jj)) = 0; 
end
for nn = 1:N_boundary
    ii = inds_pipe_boundary(nn,1)+1;
    jj = inds_pipe_boundary(nn,2)+1;
    perimiter = boundary_perimiter_span(nn);
    N_faces = VS_boundaries(ii-1,jj-1); % Number of faces the cell has which comprise the interface
    h_star = h*perimiter/N_faces; % modified heat transfer coefficient based on the ratio of the "actual" surface area of the cylidner to the surface area of the discretized interface.
    
    % At cell location
    K1 = -k/dxy;
    q = q+1;
    iii(q) = G(ii,jj);  jjj(q) = G(ii,jj);  s(q) = K1; % At cell
    
    % check left boundary
    K3 = 0;
    if ~VS_pipe(ii-2,jj-1)
        ro = ((ii-2-0.5-x_center).^2+(jj-1-0.5-y_center).^2)^0.5;
        % These values obtained from reformulating a boundary condition of the third kind into a*dT/dx+b*T=c
        K2 = k/dxy - h_star*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))); % Note regular h in denominator, not h_star!
        K3 = K3 + (-h_star*T_b*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))))/N_faces; % K3 is an average of K3 obtained for each interface
    else
        K2 = 0;
    end
    q=q+1;
    iii(q) = G(ii,jj);  jjj(q) = G(ii-1,jj);  s(q) = K2/N_faces; % left of pipe
    % check right boundary
    if ~VS_pipe(ii,jj-1)
        ro = ((ii-0.5-x_center).^2+(jj-1-0.5-y_center).^2)^0.5;
        K2 = k/dxy - h_star*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))); % Note regular h in denominator, not h_star!
        K3 = K3 + (-h_star*T_b*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))))/N_faces; % K3 is an average of K3 obtained for each interface
    else
        K2 = 0;
    end
    q=q+1;
    iii(q) = G(ii,jj); jjj(q) = G(ii+1,jj); s(q) = K2/N_faces; % right of pipe
    % check bottom boundary
    if ~VS_pipe(ii-1,jj-2)
        ro = ((ii-1-0.5-x_center).^2+(jj-2-0.5-y_center).^2)^0.5;
        K2 = k/dxy - h_star*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))); % Note regular h in denominator, not h_star!
        K3 = K3 + (-h_star*T_b*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))))/N_faces; % K3 is an average of K3 obtained for each interface
    else
        K2 = 0;
    end
        q = q+1;
        iii(q) = G(ii,jj);  jjj(q) = G(ii,jj-1);  s(q) = K2/N_faces; % bottom of pipe
    
    % check top boundary
    if ~VS_pipe(ii-1,jj)
        ro = ((ii-1-0.5-x_center).^2+(jj-0.5-y_center).^2)^0.5;
        K2 = k/dxy - h_star*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))); % Note regular h in denominator, not h_star!
        K3 = K3 + (-h_star*T_b*(1-log(ro/ri)/(k./h/dxy/ri+log(ro/ri))))/N_faces; % K3 is an average of K3 obtained for each interface  
    else
        K2 = 0;
    end
    q=q+1;
    iii(q) = G(ii,jj); jjj(q) = G(ii,jj+1); s(q) = K2/N_faces; % top of pipe
    RHS_BC_pipe(G(ii,jj)) = K3;
    
end
% Build the sparse matrix of the boundary conditions
M_pipe = sparse(iii(1:q), jjj(1:q), s(1:q), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));


end

