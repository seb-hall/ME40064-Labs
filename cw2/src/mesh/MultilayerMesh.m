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
        layers    MeshLayer % array of layer properties
        total_density       double
    end

    methods

        %% Mesh constructor
        function obj = MultilayerMesh(xmin, xmax, element_count, order, D, lambda, layers)
            
            obj = obj@Mesh(xmin, xmax, element_count, order, D, lambda);
            obj.layers = layers;


            % recalculate nodes and element counts based on layer densities
            obj.total_density = 0.0;

            for l = 1:length(layers)
                obj.total_density = obj.total_density + layers(l).density_ratio;
            end

            obj.element_count = 0;

            for l = 1:length(layers)

                layer_density = obj.layers(l).density_ratio;
                layer_element_count = round((layer_density / obj.total_density) * element_count);

                obj.layers(l).element_count = layer_element_count;
                obj.layers(l).layer_offset = obj.element_count + 1; % starting index for this layer

                obj.element_count = obj.element_count + layer_element_count;
            end

            obj.node_count = (obj.element_count * order) + 1;
            obj.node_coords = zeros(1, obj.node_count);

            obj.elements = MeshElement.empty(obj.element_count, 0);
            
        end

        function obj = Generate(obj)

            disp('Generating multilayer mesh...');

            % generate per-layer uniform node coordinates

            current_node = 1;

            for l = 1:length(obj.layers)
                
                % calculate layer xmin and xmax

                layer_xmin = obj.layers(l).x;

                if l < length(obj.layers)
                    layer_xmax = obj.layers(l + 1).x;
                else
                    layer_xmax = obj.xmax;
                end 

                layer_nodes = obj.layers(l).element_count * obj.order + 1;

                layer_coords = linspace(layer_xmin, layer_xmax, layer_nodes);


                if l == 1
                    nodes_to_add = layer_coords;
                else
                    nodes_to_add = layer_coords(2:end);  % Skip duplicate boundary node
                end

                % Add nodes to global coordinate array
                num_new_nodes = length(nodes_to_add);
                obj.node_coords(current_node : current_node + num_new_nodes - 1) = nodes_to_add;
                current_node = current_node + num_new_nodes;

                
            end

            % Generate elements
            for e = 1:obj.element_count

                % Determine global node IDs for this element
                node_start = (e - 1) * obj.order + 1;
                node_ids = node_start:(node_start + obj.order);

                % Coordinates for this element
                coords = obj.node_coords(node_ids);

                midpoint = (coords(1) + coords(end)) / 2;

                % Determine which layer this element is in
                layer_index = 1;
                for l = 1:length(obj.layers)
                    if midpoint >= obj.layers(l).x
                        layer_index = l;
                    end
                end

                D = obj.layers(layer_index).D;
                lambda = -(obj.layers(layer_index).beta + obj.layers(layer_index).gamma);

                % Create MeshElement object
                obj.elements(e) = MeshElement(node_ids, coords, obj.order, D, lambda);
            end
        end

    end
end

