function solution = DiffusionReactionSolver(mesh, D, lambda, leftBoundary, rightBoundary)
%% DiffusionReactionSolver - solves the steady-state diffusion-reaction equation for a 1D mesh
% % Inputs:
%   mesh - 1D finite element mesh structure
%   D - diffusion coefficient
%   lambda - reaction rate
%   leftBoundary - left boundary condition (BoundaryCondition object)
%   rightBoundary - right boundary condition (BoundaryCondition object)
%
% % Outputs:
%   solution - solution vector at mesh nodes

    % calculate number of nodes and elements
    Ne = mesh.ne;
    Nn = Ne + 1;

    % initialise global matrix
    globalMatrix = zeros(Nn, Nn);

    % assemble global matrix
    for eID = 1:Ne

        % get element matrices - reaction is redundant for laplace but included for completeness
        diffusionElementMatrix = DiffusionElemMatrix(D, eID, mesh);
        reactionElementMatrix = ReactionElemMatrix(lambda, eID, mesh);

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

    % here we would assemble the source vector if there were any source terms
    % however, for now we assume there are none, so it remains zero

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
end

