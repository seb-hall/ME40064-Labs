classdef ElementMatrices

    methods (Static)

        function matrix = DiffusionElemMatrix(D, elem_size)
            
            % create base matrix
            matrix = [1, -1; -1, 1];    

            % apply matrix scaling
            matrix = matrix * (D / elem_size);
        
        end

        function matrix = ReactionElemMatrix(lambda, elem_size)

            % create base matrix
            matrix = [2, 1; 1, 2];    

            % apply matrix scaling
            matrix = matrix * (lambda * elem_size / 6);
        
        end

        function matrix = MassElemMatrix(elem_size)

            % create base matrix
            matrix = [2, 1; 1, 2];    

            % apply matrix scaling
            matrix = matrix * (elem_size / 6);
        
        end
    end
end