%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07

classdef VoxelSpace < matlab.mixin.Copyable
    %VOXELSPACE Stores all physical info about a voxel space. 
    % Note that this is a copyable handle object. Calling copy(VoxelSpace) creates a shallow copy.
    %
    properties
        opaque_voxels           % (3D logical):                 True if voxels are opaque
        opaque_emissivities     % (3D double) [-]:              Emissivities of opaque voxels
        surface_normals         % (3D cell):                    1x3 surface normal for each opaque surface voxel
        surface_areas           % (3D double) [vx^2]:           Surface area of each opaque surface voxel
        PM_voxels               % (3D logical):                 True for participating media voxels
        PM_absorption_coeffs    % (3D double) [1/vx]:           Linear absorption coefficients of PM voxels
        refractive_indexes      % (3D double):                  Refractive index of voxel or of the medium that a surface voxel emits into!
        thermal_conductivity    % (3D double/scalar) [W/(m-K)]: Thermal conductivity of each voxel (only used in conduction problem)
        density                 % (3D double/scalar) [kg/m^3]:  Density of each voxel (only used in transient problem)
        specific_heat           % (3D double) [J/(kg-K)]:       Specific heat of each voxel (only used in transient problem)
        reflective_BCs          % (2x3 logical):                Defines which boundaries of the voxel space are specularly reflective
                                %                               Rows are lower/upper bound and columns are X,Y,Z axis
        size                    % (1x3 double (int)):           Size in each dimension (X,Y,Z) of voxel space
        voxel_scale             % (scalar double) [m/vx]:       Physical scale of voxel space.
        fluxes                  % (1D Flux):                    1D vector of Flux objects
        wavelength_band         % (1x2 or empty [um]):          If empty, assume properties are for entire spectrum (a gray voxel space)
        ns_normals;             % (scalar double (int)):        Neighborhood size for normal calculation
        end
    methods
        % Constructor
        function obj = VoxelSpace()
        end

        function addFlux(obj,flux)
            if ~isa(flux,"Flux")
                error("Must be flux object")
            end
            obj.fluxes{end+1} = flux;
        end
    end
end

