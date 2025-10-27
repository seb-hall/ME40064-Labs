function solution = UniformLaplaceSolver(xmin, xmax, Ne)
    
    % create 1D uniform mesh
    mesh = OneDimLinearMeshGen(xmin, xmax, Ne);

    % set parameters for Laplace equation
    lambda = 0;
    D = 1; 

    % calculate number of nodes
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
    
    % assemble source vector - for a uniform mesh we could do this in a more simple way
    % but this will work for non-uniform meshes too
    for eID = 1:Ne

        % get element size
        elemSize = mesh.elem(eID).x(2) - mesh.elem(eID).x(1);
        
        % insert into source vector
        sourceVector(eID) = sourceVector(eID) + elemSize / 2;
        sourceVector(eID + 1) = sourceVector(eID + 1) + elemSize / 2;
    end

    % Apply Neumann boundary condition
    sourceVector(1) =  sourceVector(1) + 2; % dc/dx(c=0) = 2 -> flux at left boundary

    % set Dirichlet boundary conditions
    globalMatrix(Nn, :) = 0; % zero out last row
    globalMatrix(:, Nn) = 0; % zero out last column
    globalMatrix(Nn, Nn) = 1; % set diagonal to 1
    sourceVector(Nn) = 0; % set value at last node to 0

    % solve system of equations
    solution = globalMatrix \ sourceVector;
    
end