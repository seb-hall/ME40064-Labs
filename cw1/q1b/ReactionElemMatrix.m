function matrix = ReactionElemMatrix(lambda, eID, msh)
%% ReactionElemMatrix - calculates a 2x2 element matrix for the linear
%% reaction operator in a 1d finite element mesh.
    
    % create base matrixs
    matrix = [2, 1; 1, 2];    

    % calulate element size
    elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

    % apply matrix scaling
    matrix = matrix * (lambda * elemSize / 6);
 
end