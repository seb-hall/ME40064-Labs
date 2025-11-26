%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  MeshElement.m
% Author       :  samh25
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A class defining a one-dimensional mesh element
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MeshElement

    properties

        order       uint8   % polynomial order (1 = linear, 2 = quadratic)
        node_ids    uint64  % global node IDs
        node_coords double  % node coordinates
        jacobian    double  % element jacobian d(x)/d(xi)
        D           double  % diffusion coefficient
        lambda      double  % reaction coefficient
    end

    methods

        %% MeshElement constructor
        function obj = MeshElement(ids, coords, order, D, lambda)
            
            % assign properties
            obj.node_ids = ids;
            obj.node_coords = coords;
            obj.order = order;
            obj.D = D;
            obj.lambda = lambda;

            % linear mapping from [-1, 1] to [x1, x2]
            % jacobian = dx/dxi = (x2 - x1) / 2
            obj.jacobian = (coords(end) - coords(1)) / 2;

        end
    end
end

