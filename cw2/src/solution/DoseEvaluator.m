%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Evaluation.m
% Author       :  samh25
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A static class defining a dose evaluator of 
%                 a solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef DoseEvaluator

    methods (Static)

        function K = EvaluateSolution(solution, target_x, c_threshold, dt)
            
            % first, find the closest node to the target
            node_index = 0;

            for i = 1:solution.mesh.node_count
                x = solution.mesh.node_coords(i);
                if x >= target_x
                    node_index = i;
                    break;
                end
            end

            fprintf("node index %d\n", node_index);
            c = solution.values(node_index, :);

            effective_t_index = 0;

            for i = 1:length(c)
                if c(i) >= c_threshold
                    effective_t_index = i;
                    break
                end
            end

            if effective_t_index == 0
                K = 0; % never exceeds threshold
                return;
            end 

            t_effective = solution.time(effective_t_index);

            % integrate concentration over time until effective_t_index
            time_range = effective_t_index:length(solution.time);
            K = trapz(solution.time(time_range), c(time_range));
        end

        function [c_dose_min, dose_vals, kappa_vals] = FindMinimumDose(mesh, tmax, dt, theta, integration_method, target_x, c_threshold, K_target)

            % Binary search for minimum effective dose - high and low starting bounds
            c_dose_low = 0;
            c_dose_high = 100; 

            tolerance = 0.1;

            dose_vals = [];
            kappa_vals = [];

            while (c_dose_high - c_dose_low) > tolerance
                c_dose_test = (c_dose_low + c_dose_high) / 2;
                
                fprintf("Testing dose = %.2f...\n", c_dose_test);
                
                % Run simulation with this dose
                lhs_boundary = BoundaryCondition();
                lhs_boundary.Type = BoundaryType.Dirichlet;
                lhs_boundary.Value = c_dose_test;
                
                rhs_boundary = BoundaryCondition();
                rhs_boundary.Type = BoundaryType.Dirichlet;
                rhs_boundary.Value = 0.0;
                
                solution = NumericSolver.SolveNumeric(...
                    mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, ...
                    @(~, ~) 0.0, integration_method);
                
                % Evaluate K at target location
                K = DoseEvaluator.EvaluateSolution(solution, target_x, ...
                    c_threshold, dt);
                
                if K > K_target
                    c_dose_high = c_dose_test;
                    fprintf("K = %.2f > %.0f: dose too high\n", K, K_target);
                else
                    c_dose_low = c_dose_test;
                    fprintf("K = %.2f < %.0f: dose too low\n", K, K_target);
                end

                dose_vals = [dose_vals, c_dose_test];
                kappa_vals = [kappa_vals, K];
            end
            
            c_dose_min = c_dose_high;  % Use upper bound to ensure K > target

        end
    end
end