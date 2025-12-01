%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Mesh.m
% Author       :  11973
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A class defining a one-dimensional mesh for
%                 finite element analysis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mesh < handle
    % inherit from handle to allow pass-by-reference

    properties

        xmin double
        xmax double
        dx double

        order double

        D double % diffusion coefficient
        lambda double % reaction coefficient

        node_count uint64
        node_coords double % coordinates of global nodes

        element_count uint64
        elements MeshElement % array of mesh elements

    end

    methods

        function obj = Mesh(xmin, xmax, element_count, order, D, lambda)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Mesh()
        %
        % Arguments:    parameters to initialise
        % Returns:      Mesh handle
        %
        % Description:  Initialises a one-dimensional mesh object
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            obj.xmin = xmin;
            obj.xmax = xmax;
            obj.dx = (xmax - xmin) / element_count;
            obj.D = D;
            obj.lambda = lambda;

            obj.order = order;

            % total number of nodes
            obj.node_count = (element_count * order) + 1;
            obj.node_coords = zeros(1, obj.node_count);

            obj.element_count = element_count;
            obj.elements = MeshElement.empty(element_count, 0);

        end

        function obj = Generate(obj)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Generate()
        %
        % Arguments:    Mesh handle
        % Returns:      Mesh handle
        %
        % Description:  Generates the mesh for the given object
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp('Generating normal mesh...');

            % generate uniform node coordinates
            obj.node_coords = linspace(obj.xmin, obj.xmax, obj.node_count);

            % generate elements
            for e = 1:obj.element_count

                % determine global node IDs for this element
                node_start = (e - 1) * obj.order + 1;
                node_ids = node_start:(node_start + obj.order);

                %  coordinates for this element
                coords = obj.node_coords(node_ids);

                % create MeshElement object
                obj.elements(e) = MeshElement(node_ids, coords, obj.order, obj.D, obj.lambda);
            end
        end

    end
end
