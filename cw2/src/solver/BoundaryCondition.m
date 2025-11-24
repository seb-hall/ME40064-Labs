classdef BoundaryCondition
    properties
        Type BoundaryType % Boundary condition type (Dirichlet or Neumann)
        Value double % Boundary condition value for Dirichlet
        ValueFunction function_handle % Boundary condition function for Neumann - parameter t, return double
    end
end