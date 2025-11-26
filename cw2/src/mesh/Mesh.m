%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Mesh.m
% Author       :  samh25
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

        order uint8
    
        node_count uint64
        node_coords double % coordinates of global nodes

        element_count uint64
        elements MeshElement % array of mesh elements

    end

    methods

        %% Mesh constructor
        function obj = Mesh(xmin, xmax, element_count, order, D, lambda)

            obj.xmin = xmin;
            obj.xmax = xmax;
            obj.element_count = element_count;
            obj.order = order;

            % total number of nodes
            obj.node_count = (element_count * order) + 1;
     
            % generate uniform node coordinates
            obj.node_coords = linspace(xmin, xmax, obj.node_count);

            % generate elements
            for e = 1:element_count

                % determine global node IDs for this element
                node_start = (e - 1) * order + 1;
                node_ids = node_start:(node_start + order);

                %  coordinates for this element
                coords = obj.node_coords(node_ids);

                % create MeshElement object
                obj.elements(e) = MeshElement(node_ids, coords, order, D, lambda);
            end
            
        end

    end
end

