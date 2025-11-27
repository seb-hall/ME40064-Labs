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
                if c(i) > c_threshold
                    effective_t_index = i;
                    break
                end
            end

            fprintf("effective t index %d\n", effective_t_index);

            if effective_t_index == 0
                K = 0; % never exceeds threshold
                return;
            end

            % integrate concentration over time until effective_t_index
            time_range = effective_t_index:length(solution.time);
            K = trapz(c(time_range)) * dt;
        end
    end
end