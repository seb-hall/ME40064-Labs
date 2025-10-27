function solution = DiffusionReactionSolver(mesh, D, lambda, leftBoundary, rightBoundary)
    
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
        elemMatrix = diffusionElementMatrix + reactionElementMatrix;
        
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
            sourceVector(1) = sourceVector(1) - leftBoundary.Value;
        case BoundaryType.Dirichlet
            globalMatrix(1, :) = 0; % zero out first row
            globalMatrix(:, 1) = 0; % zero out first column
            globalMatrix(1, 1) = 1; % set diagonal to 1
            sourceVector(1) = leftBoundary.Value; % set value at first node            
    end

    % Apply right boundary condition
    switch rightBoundary.Type
        case BoundaryType.Neumann
            sourceVector(Nn) = sourceVector(Nn) - rightBoundary.Value;
        case BoundaryType.Dirichlet
            globalMatrix(Nn, :) = 0; % zero out last row
            globalMatrix(:, Nn) = 0; % zero out last column
            globalMatrix(Nn, Nn) = 1; % set diagonal to 1
            sourceVector(Nn) = rightBoundary.Value; % set value at last node            
    end

    % solve system of equations
    solution = globalMatrix \ sourceVector;
end

