%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  NumericSolverTest.m
% Author       :  11973
% Created      :  2025-11-27 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Test suite for NumericSolver class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = NumericSolverTest
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
        mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @SourceFunction, integration_method);

    % analytical solution
    t_analytic = 0:dt:tmax;
    c_exact = 10 * (1 - exp(lambda * t_analytic));

    % chose random node (3 in this case) to compare
    c_numeric = numeric_solution.values(3,:);
    error = norm(c_numeric - c_exact) / norm(c_exact);

    tolerance = 1e-3;
    verifyLessThan(testCase, error, tolerance);
end

function s = SourceFunction(x, t)
    s = 10;
end
