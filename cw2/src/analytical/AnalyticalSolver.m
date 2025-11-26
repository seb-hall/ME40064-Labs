%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  AnalyticalSolver.m
% Author       :  samh25
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A static class defining an analytical solver
%                 for the transient diffusion equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef AnalyticalSolver

    methods (Static)

        function solution = SolveAnalytical(mesh, tmax, dt)
            
            % time vector
            time_vector = 0:dt:tmax;
            solution = Solution(mesh, time_vector);

            % loop over time steps
            for step = 1:length(time_vector)

                t = time_vector(step);
                timestep_results = zeros(1, mesh.node_count);
                
                % loop over nodes
                for i = 1:mesh.node_count
                    x = mesh.node_coords(i);
                    timestep_results(i) = TransientAnalyticSoln(x, t);
                end

                solution.set_values(timestep_results, step);
            end

        end
    end
end