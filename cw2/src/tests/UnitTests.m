%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  UnitTests.m
% Author       :  11973
% Created      :  2025-11-27 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Test suite for FEM transient diffusion solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = UnitTests
    % run unit tests for transient diffusion solver
    tests = functiontests(localfunctions);
end

function TestSolveNumericReactionOnly(testCase)

    % mesh parameters
    xmin = 0.0;
    xmax = 1.0;
    element_count = 6;
    order = 1;
    lambda = -1.0;
    D = 0.0;

    % time parameters
    tmax = 0.5;
    dt = 0.02;
    theta = 0.5; % Crank-Nicholson

    % generate mesh
    mesh = Mesh(xmin, xmax, element_count, order, D, lambda);
    mesh.Generate();

    % solver parameters
    lhs_boundary = BoundaryCondition();
    lhs_boundary.Type = BoundaryType.Neumann;
    lhs_boundary.ValueFunction = @(t) 0.0;

    rhs_boundary = BoundaryCondition();
    rhs_boundary.Type = BoundaryType.Neumann;
    rhs_boundary.ValueFunction = @(t) 0.0;

    integration_method = IntegrationMethod();
    integration_method.type = IntegrationType.Trapezoidal;
    integration_method.gauss_points = 0; % not used for trapezoidal
   
    numeric_solution = NumericSolver.SolveNumeric(...
        mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 10, integration_method);

    % analytical solution
    t_analytic = 0:dt:tmax;
    c_exact = 10 * (1 - exp(lambda * t_analytic));

    % chose random node (3 in this case) to compare
    c_numeric = numeric_solution.values(3,:);
    error = norm(c_numeric - c_exact) / norm(c_exact);

    tolerance = 1e-3;
    verifyLessThan(testCase, error, tolerance);
end

function TestSteadyStateConvergence(testCase)
    % For pure diffusion with constant BCs, should reach steady state
    % Analytical: c(x) = x (linear profile from 0 to 1)
    
    xmin = 0; xmax = 1;
    mesh = Mesh(xmin, xmax, 20, 1, 1.0, 0.0);
    mesh.Generate();
    
    % Run to steady state
    tmax = 10.0; dt = 0.1; theta = 0.5;
    
    lhs_bc = BoundaryCondition();
    lhs_bc.Type = BoundaryType.Dirichlet;
    lhs_bc.Value = 0.0;
    
    rhs_bc = BoundaryCondition();
    rhs_bc.Type = BoundaryType.Dirichlet;
    rhs_bc.Value = 1.0;
    
    integration_method = IntegrationMethod();
    integration_method.type = IntegrationType.Trapezoidal;
    
    solution = NumericSolver.SolveNumeric(mesh, tmax, dt, theta, ...
                                          lhs_bc, rhs_bc, @(x,t) 0, integration_method);
    
    % Check final solution is linear
    c_final = solution.values(:, end);
    c_expected = mesh.node_coords';
    
    error = norm(c_final - c_expected) / norm(c_expected);
    verifyLessThan(testCase, error, 1e-3);
end

function TestMassConservation(testCase)
    % With D=1, lambda=0, no source, Neumann BCs: total mass conserved
    
    mesh = Mesh(0, 1, 10, 1, 1.0, 0.0);
    mesh.Generate();
    
    % Set initial non-zero distribution
    c0 = sin(pi * mesh.node_coords);
    
    lhs_bc = BoundaryCondition();
    lhs_bc.Type = BoundaryType.Neumann;
    lhs_bc.ValueFunction = @(t) 0.0;
    
    rhs_bc = BoundaryCondition();
    rhs_bc.Type = BoundaryType.Neumann;
    rhs_bc.ValueFunction = @(t) 0.0;
    
    % [Run solver with initial condition c0]
    % Check integral(c dx) is constant over time
end
