function [mesh] = OneDimSimpleRefinedMeshGen(xmin,xmax,Ne)
%%This function generates a one dimensional, equispaced, linear finite
%%element mesh, with Ne number of elements, between the points at x
%%position xmin and xmax.

    mesh.ne = Ne; %set number of elements
    mesh.ngn = Ne+1; %set number of global nodes
    mesh.nvec = zeros(mesh.ngn,1); %allocate vector to store global node values
    
    %dx = (xmax - xmin)/Ne; %calculate element size
    mesh.nvec(1) = xmin;
    mesh.nvec(mesh.ngn) = xmax;
    
    if(mesh.ngn > 1)
        for j=2:(mesh.ngn-1)
          mesh.nvec(j) = (xmax - xmin)/(2^(j-1)) + mesh.nvec(j-1)
        end
    end
    %mesh.nvec = xmin:dx:xmax;
    
    %loop over elements & set the element properties
    for i=1:Ne

        %set global IDs of the nodes
        mesh.elem(i).n(1) = i;
        mesh.elem(i).n(2) = i+1;
        
        %set spatial positions of nodes
        mesh.elem(i).x(1) = mesh.nvec(mesh.elem(i).n(1));
        mesh.elem(i).x(2) = mesh.nvec(mesh.elem(i).n(2));
        %mesh.elem(i).x(1) = xmin + (i-1)*dx;
        %mesh.elem(i).x(2) = xmin + i*dx ;

        

        %set element Jacobian based on mapping to standard element
        mesh.elem(i).J = 0.5*(mesh.elem(i).x(2) - mesh.elem(i).x(1)); %this is assuming standard element of -1 to 1
       

    end

end