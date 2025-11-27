classdef ElementMatrices

    methods (Static)

        function matrix = DiffusionElemMatrix(element, method)
            
            elem_size = element.node_coords(end) - element.node_coords(1); 

            if method.type == IntegrationType.Trapezoidal
                
                % create base matrix
                matrix = eye(element.order + 1);
                
                for i = 1:(element.order + 1)
                    for j = 1:(element.order + 1)
                        if i ~= j
                            matrix(i, j) = -1;
                        end
                    end
                end

                % apply matrix scaling
                matrix = matrix * (element.D / elem_size);
            
            else

                matrix = zeros(element.order + 1);
                [xi, wi] = ElementMatrices.GaussQuadraturePoints(method.gauss_points);

                for i = 1:length(xi)
                    dN_dxi = ElementMatrices.ShapeFunctionDerivatives(element.order, xi(i));

                    J = element.jacobian;
                    dN_dx = dN_dxi / J;

                    % compute contribution to stiffness matrix
                    matrix = matrix + (element.D * (dN_dx' * dN_dx)) * (wi(i) * J);
                end

            end
        
        end

        function matrix = ReactionElemMatrix(element, method)

            elem_size = element.node_coords(end) - element.node_coords(1); 

            if method.type == IntegrationType.Trapezoidal

                % create base matrix
                matrix = eye(element.order + 1) * 2;

                for i = 1:(element.order + 1)
                    for j = 1:(element.order + 1)
                        if i ~= j
                            matrix(i, j) = 1;
                        end
                    end
                end

                % apply matrix scaling
                matrix = matrix * (element.lambda * elem_size / 6);

            else

                matrix = zeros(element.order + 1);
                [xi, wi] = ElementMatrices.GaussQuadraturePoints(method.gauss_points);

                for i = 1:length(xi)
                    N = ElementMatrices.ShapeFunctions(element.order, xi(i));

                    J = element.jacobian;

                    % compute contribution to stiffness matrix
                    matrix = matrix + (element.lambda * (N' * N)) * (wi(i) * J);
                end
            end
        end

        function matrix = MassElemMatrix(element, method)

            elem_size = element.node_coords(end) - element.node_coords(1); 

            if method.type == IntegrationType.Trapezoidal

                % create base matrix
                matrix = eye(element.order + 1) * 2;

                for i = 1:(element.order + 1)
                    for j = 1:(element.order + 1)
                        if i ~= j
                            matrix(i, j) = 1;
                        end
                    end
                end

                % apply matrix scaling
                matrix = matrix * (elem_size / 6);
                
            else

                matrix = zeros(element.order + 1);
                [xi, wi] = ElementMatrices.GaussQuadraturePoints(method.gauss_points);

                for i = 1:length(xi)
                    N = ElementMatrices.ShapeFunctions(element.order, xi(i));

                    J = element.jacobian;

                    % compute contribution to stiffness matrix
                    matrix = matrix + (N' * N) * (wi(i) * J);
                end

            end
        
        end

        function matrix = ForceMatrix(element, method)

            elem_size = element.node_coords(end) - element.node_coords(1); 

            if method.type == IntegrationType.Trapezoidal

                % create base matrix
                matrix = ones(element.order + 1, 1);

                % apply matrix scaling
                matrix = matrix * (elem_size / 2);
                
            else

                matrix = zeros(element.order + 1, 1);
                [xi, wi] = ElementMatrices.GaussQuadraturePoints(method.gauss_points);

                for i = 1:length(xi)
                    N = ElementMatrices.ShapeFunctions(element.order, xi(i));

                    J = element.jacobian;

                    % compute contribution to stiffness matrix
                    matrix = matrix + N' * (wi(i) * J);
                end

            end

        end
    end

    methods (Static, Access = private)

        function [xi, wi] = GaussQuadraturePoints(n)

            switch n
                case 1
                    xi = 0;
                    wi = 2;
                case 2
                    xi = [-1/sqrt(3), 1/sqrt(3)];
                    wi = [1, 1];
                case 3
                    xi = [-sqrt(3/5), 0, sqrt(3/5)];
                    wi = [5/9, 8/9, 5/9];
                otherwise
                    error('Gauss quadrature for n > 3 not implemented.');
            end

        end

        function N = ShapeFunctions(order, xi)

            switch order
                case 1 % linear
                    N = [(1 - xi) / 2, (1 + xi) / 2];
                case 2 % quadratic
                    N = [xi * (xi - 1) / 2, (1 - xi^2), xi * (xi + 1) / 2];
                otherwise
                    error('Shape functions for order > 2 not implemented.');
            end

        end

        function dN_dxi = ShapeFunctionDerivatives(order, xi)

            switch order
                case 1 % linear
                    dN_dxi = [-0.5, 0.5];
                case 2 % quadratic
                    dN_dxi = [xi - 0.5, -2 * xi, xi + 0.5];
                otherwise
                    error('Shape function derivatives for order > 2 not implemented.');
            end

        end
    end
end