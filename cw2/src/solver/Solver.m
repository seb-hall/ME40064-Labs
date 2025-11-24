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

        mesh % mesh object  

        % time variant properties, all should be a function handle with time as input
        D function_handle % diffusion coefficient 
        lambda function_handle % reaction rate
        K function_handle % global stiffness matrix
        M function_handle % global mass matrix



    end

    methods

        % generic constructor
        function obj = Solver(dt, mesh, D, lambda, K, M)

            obj.t = 0;
            obj.dt = dt;
            obj.mesh = mesh;
            obj.D = D;
            obj.lambda = lambda;
            obj.K = K;
            obj.M = M;

        end

        function obj = Update(obj)

             % calculate number of nodes and elements
            Ne = obj.mesh.ne;
            Nn = Ne + 1;

            % initialise global matrix
            globalMatrix = zeros(Nn, Nn);

            % assemble global matrix
            for eID = 1:Ne

                % get element matrices - reaction is redundant for laplace but included for completeness
                diffusionElementMatrix = DiffusionElemMatrix(obj.D(obj.t), eID, obj.mesh);
                reactionElementMatrix = ReactionElemMatrix(obj.lambda(obj.t), eID, obj.mesh);

                % combine element matrices
                elemMatrix = diffusionElementMatrix - reactionElementMatrix;
                
                % insert into global matrix
                globalMatrix(eID, eID) = globalMatrix(eID, eID) + elemMatrix(1, 1);
                globalMatrix(eID, eID + 1) = globalMatrix(eID, eID + 1) + elemMatrix(1, 2);
                globalMatrix(eID + 1, eID) = globalMatrix(eID + 1, eID) + elemMatrix(2, 1);
                globalMatrix(eID + 1, eID + 1) = globalMatrix(eID + 1, eID + 1) + elemMatrix(2, 2);

            end

            % initialise source vector
            sourceVector = zeros(Nn, 1);

            % Apply left boundary condition
            switch leftBoundary.Type

                case BoundaryType.Neumann

                    % directly modify source vector for Neumann
                    sourceVector(1) = sourceVector(1) - leftBoundary.Value;

                case BoundaryType.Dirichlet

                    % apply Dirichlet condition to source vector
                    for j = 2:Nn
                        sourceVector(j) = sourceVector(j) - globalMatrix(j, 1) * leftBoundary.Value;
                    end

                    % modify global matrix
                    globalMatrix(1, :) = 0;
                    globalMatrix(:, 1) = 0;
                    globalMatrix(1, 1) = 1;

                    % set value at first node
                    sourceVector(1) = leftBoundary.Value; 
            end

            % Apply right boundary condition
            switch rightBoundary.Type

                case BoundaryType.Neumann

                    % directly modify source vector for Neumann
                    sourceVector(Nn) = sourceVector(Nn) - rightBoundary.Value;

                case BoundaryType.Dirichlet

                    % apply Dirichlet condition to source vector
                    for j = 2:(Nn-1)
                        sourceVector(j) = sourceVector(j) - globalMatrix(j, Nn) * rightBoundary.Value;
                    end

                    % modify global matrix
                    globalMatrix(Nn, :) = 0;
                    globalMatrix(:, Nn) = 0;
                    globalMatrix(Nn, Nn) = 1;

                    % set value at last node
                    sourceVector(Nn) = rightBoundary.Value; 
            end

            % solve system of equations
            solution = globalMatrix \ sourceVector;

            % update time
            obj.t = obj.t + obj.dt;
        end

        function obj = PlotResults(obj)
            % placeholder for plotting results
        end

        function obj = Test(obj)
            % placeholder for testing solver
        end
    end

end
