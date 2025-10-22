function matrix = LaplaceElemMatrix(D, eID, msh)
%% LaplaceElemMatrix - calculates a 2x2 element matrix for the diffusion 
%% operator in a 1d finite element mesh.
    
    % create base matrix
    matrix = [1, -1; -1, 1];    

    % calulate element size
    elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

    % scale base matrix by element size
    matrix = matrix * (D / elemSize);
 
end