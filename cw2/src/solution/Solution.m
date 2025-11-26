%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Solution.m
% Author       :  samh25
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A class defining a solution to the transient
%                 diffusion equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Solution < handle
    % inherit from handle to allow pass-by-reference behaviour

    properties

        mesh        Mesh    % handle to mesh object
        time        double  % time series - 1 x Nsteps
        values      double  % solution values - Nnodes x Nsteps

    end

    methods

        %% Solution constructor
        function obj = Solution(mesh, time_vector)
            
            % assign properties
            obj.mesh = mesh;
            obj.time = time_vector;
            obj.values = zeros(mesh.node_count, length(time_vector));
            
        end

        %% Set solution values at given time step
        function SetValues(obj, values, step)
            % set solution values at given time step
            obj.values(:, step) = values(:);
        end
    end
end

