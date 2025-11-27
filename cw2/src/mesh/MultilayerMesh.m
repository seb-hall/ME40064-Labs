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


classdef MultilayerMesh < Mesh
    % inherit from handle to allow pass-by-reference

    properties
        layer_properties LayerProperties % array of layer properties
    end

    methods

        %% Mesh constructor
        function obj = MultilayerMesh(xmin, xmax, element_count, order, D, lambda, layer_properties)
            
            obj = obj@Mesh(xmin, xmax, element_count, order, D, lambda);
            obj.layer_properties = layer_properties;
            
        end

        function obj = Generate(obj)

            disp('Generating multilayer mesh...');

            % generate uniform node coordinates
            obj.node_coords = linspace(obj.xmin, obj.xmax, obj.node_count);

            % generate elements
            for e = 1:obj.element_count

                % determine global node IDs for this element
                node_start = (e - 1) * obj.order + 1;
                node_ids = node_start:(node_start + obj.order);

                %  coordinates for this element
                coords = obj.node_coords(node_ids);

                midpoint = (coords(1) + coords(end)) / 2;

                % determine which layer this element is in
                layer_index = 1;

                for l = 1:length(obj.layer_properties)
                    if midpoint >= obj.layer_properties(l).x
                        layer_index = l;
                    end
                end

                D = obj.layer_properties(layer_index).D;
                lambda = -(obj.layer_properties(layer_index).beta + obj.layer_properties(layer_index).gamma);

                % create MeshElement object
                obj.elements(e) = MeshElement(node_ids, coords, obj.order, D, lambda);
            end
        end

    end
end

