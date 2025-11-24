%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Solver.m
% Author       :  samh25
% Created      :  2025-11-24 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Class definition for generic solver for the 
%                 transient diffusion-reaction equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef Solver < handle

    properties 
        t double % current time
        dt double % time step size
        
        step_count uint64 % current step count
        steps uint64 % total steps
        
        % weighting factor for time integration
        % 0 = Forward Euler, 0.5 = Crank-Nicholson, 1 = Backward Euler
        theta double 

        mesh % mesh object  
        solution % global solution vector

        % physical properties
        D double
        lambda double

        % global matrices
        M double
        K double

        % boundary conditions
        left_boundary BoundaryCondition
        right_boundary BoundaryCondition

        % source function - take x and t, return value (double)
        source_function function_handle

    end

    methods

        % generic constructor
        function obj = Solver(dt, steps, theta, mesh, D, lambda, left_boundary, right_boundary, source_function)

            obj.t = 0;
            obj.dt = dt;
            obj.steps = steps;
            obj.step_count = 1; % MATLAB indexing starts at 1

            obj.theta = theta;

            obj.mesh = mesh;
            obj.solution = zeros(mesh.ngn, obj.steps + 1);

            obj.D = D;
            obj.lambda = lambda;

            obj.left_boundary = left_boundary;
            obj.right_boundary = right_boundary;

            obj.source_function = source_function;
            
            [obj.K, obj.M] = obj.CreateStiffnessMatrix();

        end

        function obj = Update(obj)

            % prevent
            if obj.step_count > obj.steps
                return;
            end
            
            c_current = obj.solution(:, obj.step_count);

            % ASSEMBLE SYSTEM MATRIX AND RHS VECTOR
            system_matrix = obj.M + obj.theta * obj.dt * obj.K;
            rhs_vector = (obj.M - (1 - obj.theta) * obj.dt * obj.K) * c_current;

            % APPLY SOURCE TERMS
            f_current = obj.CreateSourceVector(obj.t);
            f_next = obj.CreateSourceVector(obj.t + obj.dt);
            rhs_vector = rhs_vector + obj.dt * (obj.theta * f_next + (1 - obj.theta) * f_current);

            % APPLY BOUNDARY CONDITIONS
            [system_matrix, rhs_vector] = obj.ApplyBoundaryConditions(system_matrix, rhs_vector);
    

            c_next = system_matrix \ rhs_vector;
            obj.solution(:, obj.step_count + 1) = c_next;

            % update time
            obj.t = obj.t + obj.dt;
            obj.step_count = obj.step_count + 1;
        end

        function [K, M] = CreateStiffnessMatrix(obj)

            Ne = obj.mesh.ne;
            Nn = obj.mesh.ngn;

            % initialise global matrix (use sparse for efficiency with large systems)
            K = sparse(Nn, Nn);
            M = sparse(Nn, Nn);

            % assemble global matrix
            for eID = 1:Ne

                % get element matrices - reaction is redundant for laplace but included for completeness
                diffM = ElementMatrices.DiffusionElemMatrix(obj.D, eID, obj.mesh);
                reactM = ElementMatrices.ReactionElemMatrix(obj.lambda, eID, obj.mesh);

                % combine element matrices
                elemK = diffM - reactM;
                
                % insert into global stiffness matrix
                K(eID, eID) = K(eID, eID) + elemK(1, 1);
                K(eID, eID + 1) = K(eID, eID + 1) + elemK(1, 2);
                K(eID + 1, eID) = K(eID + 1, eID) + elemK(2, 1);
                K(eID + 1, eID + 1) = K(eID + 1, eID + 1) + elemK(2, 2);
                
                % get mass element matrix
                elemM = ElementMatrices.MassElemMatrix(obj.D, eID, obj.mesh);

                % insert into global mass matrix
                M(eID, eID) = M(eID, eID) + elemM(1, 1);
                M(eID, eID + 1) = M(eID, eID + 1) + elemM(1, 2);
                M(eID + 1, eID) = M(eID + 1, eID) + elemM(2, 1);
                M(eID + 1, eID + 1) = M(eID + 1, eID + 1) + elemM(2, 2);

            end
        end 

        function [lhs, rhs] = ApplyBoundaryConditions(obj, lhs, rhs)

            % apply left boundary condition

            switch obj.left_boundary.Type

                case BoundaryType.Dirichlet
                    lhs(1, :) = 0; % clear row
                    lhs(1, 1) = 1; % set diagonal to 1
                    rhs(1) = obj.left_boundary.Value; % set value

                case BoundaryType.Neumann
                    rhs(1) = rhs(1) + obj.left_boundary.ValueFunction(obj.t); % apply flux
            
            end

            switch obj.right_boundary.Type
                 case BoundaryType.Dirichlet
                    lhs(end, :) = 0; % clear row
                    lhs(end, end) = 1; % set diagonal to 1
                    rhs(end) = obj.right_boundary.Value; % set value

                case BoundaryType.Neumann
                     rhs(end) = rhs(end) + obj.right_boundary.ValueFunction(obj.t); % apply flux
            end

        end

        function F = CreateSourceVector(obj, t)

            F = zeros(obj.mesh.ngn, 1);

            % return if no source function defined
            if (~ishandle(obj.source_function))
                return;
            end

            for i = 1:obj.mesh.ne
                
                h = obj.mesh.elem(i).x(2) - obj.mesh.elem(i).x(1);

                midpoint = (obj.mesh.elem(i).x(1) + obj.mesh.elem(i).x(2)) / 2;
                f_val = obj.source_function(midpoint, t);

                % Local Force Vector for linear element (Int N^T * s dx)
                f_local = f_val * h / 2 * [1; 1];

                nodes = [i, i + 1];
                F(nodes) = F(nodes) + f_local;
            
            end




        end
    
    end

end
