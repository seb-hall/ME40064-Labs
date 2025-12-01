%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  BoundaryCondition.m
% Author       :  11973
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A class defining a boundary condition for the
%                 transient diffusion reaction equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BoundaryCondition

    properties

        Type BoundaryType % Boundary condition type (Dirichlet or Neumann)
        Value double % Boundary condition value for Dirichlet
        ValueFunction function_handle % Boundary condition function for Neumann - parameter t, return double
        
    end

end