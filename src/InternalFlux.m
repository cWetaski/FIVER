%   AUTHOR: Charles Wetaski
%   LAST CHECKED: 2024-06-07

classdef InternalFlux < Flux
    %INTERNALFLUX Defines an internal flux within a voxel space
    %   Generates rays according to a position and direction distribution of a source flux within a voxel space as 
    %   well as the overallpower (in W) associated with the flux
    %
    % Properties inherited from Flux.m
    methods
        function obj = InternalFlux(power,ray_gen_function)
            %EXTERNALFLUX Constructs an ExternalFlux Object and assigns it to the VoxelSpace
            obj.power = power;
            obj.ray_gen_function = ray_gen_function;
        end
        
        function [rays_pos,rays_dir] = GenerateRays(obj,N_rays)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            rays = obj.ray_gen_function(N_rays);
            rays_pos = rays(:,1:3);
            rays_dir = rays(:,4:6);
        end
    end
end

