function [mesh] = OneDimMeshGen(xmin, xmax, Ne, order)


    mesh.ne = Ne; %set number of elements
    mesh.order = order;

    mesh.ngn = Ne * (order - 1) + 1;

    mesh.nvec = linspace(xmin, xmax, mesh.ngn);

    elementLength = (xmax - xmin) / Ne;
    stride = order - 1;

    %loop over elements & set the element properties
    for i=1:Ne

        startNodeIndex = (i - 1) * stride + 1;
        nodeIndices = startNodeIndex : (startNodeIndex + stride);

        mesh.elem(i).n = nodeIndices;
        mesh.elem(i).x = mesh.nvec(nodeIndices);

        mesh.elem(i).J = 0.5 * elementLength;
    end

end
