%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  L2Error.m
% Author       :  samh25
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A class defining a solution to the transient
%                 diffusion equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef L2Error < handle
    % inherit from handle to allow pass-by-reference behaviour

    properties

        ref_solution    Solution % handle to reference solution object
        num_solution    Solution % handle to solution object

        time            double   % time series - 1 x Nsteps
        l2_error        double   % L2 error at each time step - 1 x Nsteps

    end

    methods

        %% Solution constructor
        function obj = L2Error(ref_solution, num_solution)
            
            % assign properties
            obj.ref_solution = ref_solution;
            obj.num_solution = num_solution;
            obj.time = ref_solution.time;

            if ref_solution.mesh.node_count ~= num_solution.mesh.node_count
                error('Reference and error solutions must have the same number of nodes');
            end

            if length(ref_solution.time) ~= length(num_solution.time)
                error('Reference and error solutions must have the same number of time steps');
            end 

            step_count = length(ref_solution.time);
            
            obj.l2_error = zeros(1, step_count);

            for step = 1:step_count
                c_ref = ref_solution.values(:, step);
                c_num = num_solution.values(:, step);
                x = ref_solution.mesh.node_coords;

                integrand = (c_ref - c_num).^2;
                obj.l2_error(step) = sqrt(trapz(x, integrand));
            end

        end

    end
end

