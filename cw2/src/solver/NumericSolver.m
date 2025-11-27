%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  NumericSolver.m
% Author       :  samh25
% Created      :  2025-11-24 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Class definition for generic solver for the 
%                 transient diffusion-reaction equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef NumericSolver

    methods (Static)

       function solution = SolveNumeric(mesh, tmax, dt, theta, left_boundary, right_boundary, source_fn, integration_method)
            
            % time vector
            time_vector = 0:dt:tmax;
            solution = Solution(mesh, time_vector);

            % === SET INITIAL CONDITION EXPLICITLY ===
            c0 = zeros(mesh.node_count, 1);
            solution.SetValues(c0, 1);  % column 1 = t=0

            [K, M] = NumericSolver.CreateGlobalMatrices(mesh, integration_method);

            % loop over time steps
            for step = 1:length(time_vector) - 1

                c_next = NumericSolver.SolveStep(mesh, solution, step, dt, theta, K, M, left_boundary, right_boundary, source_fn, integration_method);
                solution.SetValues(c_next, step + 1);

            end

        end

        function c = SolveStep(mesh, solution, step, dt, theta, K, M, left_boundary, right_boundary, source_fn, integration_method)

            c_current = solution.values(:, step);

            t = (step - 1) * dt; % current time, converted to 0-based index

            % assemble system matrix and rhs vector
            system_matrix = M + theta * dt * K;
            rhs_vector = (M - (1 - theta) * dt * K) * c_current;

            % add source term
            f_current = NumericSolver.CreateSourceVector(mesh, t, source_fn, integration_method);
            f_next = NumericSolver.CreateSourceVector(mesh, t + dt, source_fn, integration_method);
            rhs_vector = rhs_vector + dt * (theta * f_next + (1 - theta) * f_current);

            % apply boundary conditions
            [system_matrix, rhs_vector] = NumericSolver.ApplyBoundaryConditions(system_matrix, rhs_vector, t + dt, left_boundary, right_boundary);

            % solve system
            c = system_matrix \ rhs_vector;

        end

        %% Create global stiffness and mass matrices
        function [K, M] = CreateGlobalMatrices(mesh, integration_method)

            num_elements = mesh.element_count;
            num_nodes = mesh.node_count;

            % initialise global matrix (use sparse for efficiency with large systems)
            K = sparse(num_nodes, num_nodes);
            M = sparse(num_nodes, num_nodes);

            for element_id = 1:num_elements

                element = mesh.elements(element_id);
                nodes = element.node_ids;
                local_size = length(nodes);

                diff_matrix = ElementMatrices.DiffusionElemMatrix(element, integration_method);
                react_matrix = ElementMatrices.ReactionElemMatrix(element, integration_method);

                k_matrix = diff_matrix - react_matrix;

                elem_size = element.node_coords(end) - element.node_coords(1); 
                m_matrix = ElementMatrices.MassElemMatrix(element, integration_method);

                % assemble into global matrices
                for i = 1:local_size
                    for j = 1:local_size
                        gi = nodes(i); gj = nodes(j);
                        K(gi, gj) = K(gi, gj) + k_matrix(i, j);
                        M(gi, gj) = M(gi, gj) + m_matrix(i, j);
                    end
                end

            end

        end

        %% Create source vector for given time
        function F = CreateSourceVector(mesh, t, source_fn, integration_method)

            F = zeros(mesh.node_count, 1);

            % return if no source function defined
            if (isempty(source_fn))
                return;
            end

            for element_id = 1:mesh.element_count
                
                element = mesh.elements(element_id);

                elem_size = element.node_coords(end) - element.node_coords(1); 
                midpoint = (element.node_coords(1) + element.node_coords(end)) / 2;

                f_val = source_fn(midpoint, t);

                % Local Force Vector for linear element (Int N^T * s dx)
                f_local = f_val * ElementMatrices.ForceMatrix(element, integration_method);

                nodes = element.node_ids;
                F(nodes) = F(nodes) + f_local;
            end

        end

        function [lhs, rhs] = ApplyBoundaryConditions(lhs, rhs, t, left_boundary, right_boundary)

            % apply left boundary condition

            switch left_boundary.Type

                case BoundaryType.Dirichlet
                    lhs(1, :) = 0;                  % clear row
                    lhs(1, 1) = 1;                  % set diagonal to 1
                    rhs(1) = left_boundary.Value;   % set value

                case BoundaryType.Neumann
                    rhs(1) = rhs(1) + left_boundary.ValueFunction(t); % apply flux
            
            end

            switch right_boundary.Type
                 case BoundaryType.Dirichlet
                    lhs(end, :) = 0;                % clear row
                    lhs(end, end) = 1;              % set diagonal to 1
                    rhs(end) = right_boundary.Value; % set value

                case BoundaryType.Neumann
                    rhs(end) = rhs(end) + right_boundary.ValueFunction(t); % apply flux
            end

        end

    end
end
