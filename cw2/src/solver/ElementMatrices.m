classdef ElementMatrices

    methods (Static)

        function matrix = ReactionElemMatrix(lambda, eID, msh)

            % create base matrixs
            matrix = [2, 1; 1, 2];    

            % calulate element size
            elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

            % apply matrix scaling
            matrix = matrix * (lambda * elemSize / 6);
        
        end

        function matrix = DiffusionElemMatrix(D, eID, msh)
            
            % create base matrix
            matrix = [1, -1; -1, 1];    

            % calulate element size
            elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

            % apply matrix scaling
            matrix = matrix * (D / elemSize);
        
        end


        function matrix = MassElemMatrix(~, eID, msh)

            % create base matrix
            matrix = [2, 1; 1, 2];    

            % calulate element size
            elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

            % apply matrix scaling
            matrix = matrix * (elemSize / 6);
        
        end
    end
end