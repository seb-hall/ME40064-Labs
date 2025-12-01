%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  ElementMatrices.m
% Author       :  11973
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A static class defining element matrix
%                 helper functions for the transient diffusion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ElementMatrices

    methods (Static)

        function matrix = DiffusionElemMatrix(element, method)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     DiffusionElemMatrix()
        %
        % Arguments:    element and integration method
        % Returns:      diffusion element matrix
        %
        % Description:  Computes the diffusion element matrix
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     ReactionElemMatrix()
        %
        % Arguments:    element and integration method
        % Returns:      reaction element matrix
        %
        % Description:  Computes the reaction element matrix
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     MassElemMatrix()
        %
        % Arguments:    element and integration method
        % Returns:      mass element matrix
        %
        % Description:  Computes the mass element matrix
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     ForceMatrix()
        %
        % Arguments:    element and integration method
        % Returns:      force element matrix
        %
        % Description:  Computes the force element matrix
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     GaussQuadraturePoints()
        %
        % Arguments:    number of points
        % Returns:      quadrature points and weights
        %
        % Description:  Looks up Gauss quadrature points and weights
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     ShapeFunctions()
        %
        % Arguments:    order and local coordinate
        % Returns:      shape function values
        %
        % Description:  Computes shape function values at given local
        %               coordinate
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     ShapeFunctionDerivatives()
        %
        % Arguments:    order and local coordinate
        % Returns:      shape function derivative values
        %
        % Description:  Computes shape function derivative values at 
        %               given local coordinate
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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