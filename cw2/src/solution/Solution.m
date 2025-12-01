%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Solution.m
% Author       :  11973
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

        function obj = Solution(mesh, time_vector)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Solution()
        %
        % Arguments:    parameters to initialise
        % Returns:      Solution handle
        %
        % Description:  Initialises a Solution object
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            % assign properties
            obj.mesh = mesh;
            obj.time = time_vector;
            obj.values = zeros(mesh.node_count, length(time_vector));
            
        end


        function obj = SetValues(obj, values, step)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     SetValues()
        %
        % Arguments:    object, values to set, time step index
        % Returns:      Solution handle
        %
        % Description:  Helper function to set solution values at a
        %               given time step
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % set solution values at given time step
            obj.values(:, step) = values(:);
        end
    end
end

