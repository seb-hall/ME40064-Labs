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
    end

    methods

        % generic constructor
        function obj = Solver(dt, steps, theta, mesh, D, lambda)

            obj.t = 0;
            obj.dt = dt;
            obj.steps = steps;
            obj.step_count = 1; % MATLAB indexing starts at 1

            obj.theta = theta;

            obj.mesh = mesh;
            obj.solution = zeros(mesh.ngn, obj.steps + 1);

            obj.D = D;
            obj.lambda = lambda;

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
            f_current = zeros(obj.mesh.ngn, 1);
            f_next = zeros(obj.mesh.ngn, 1);
            %rhs_vector = rhs_vector + obj.dt * (obj.theta * f_next + (1 - obj.theta) * f_current);

            % APPLY BOUNDARY CONDITIONS
            % TBD

            % --- Left Boundary (x=0) ---
            % Condition: c(0, t) = 0
            % Node Index: 1
            system_matrix(1, :) = 0;   % Clear the entire row
            system_matrix(1, 1) = 1;   % Set diagonal to 1
            rhs_vector(1)       = 0;   % Set value to 0
            
            % --- Right Boundary (x=1) ---
            % Condition: c(1, t) = 1
            % Node Index: Last Node (obj.mesh.ngn)
            lastNode = obj.mesh.ngn;
            system_matrix(lastNode, :)        = 0;   % Clear the entire row
            system_matrix(lastNode, lastNode) = 1;   % Set diagonal to 1
            rhs_vector(lastNode)              = 1;   % Set value to 1 (The driving force!)

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
    
    end

end
