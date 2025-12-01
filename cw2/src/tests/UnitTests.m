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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function:     TestSolveNumericReactionOnly()
%
% Arguments:    test case
% Returns:      none
%
% Description:  Tests NumericSolver for pure reaction case
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function:     TestSteadyStateConvergence()
%
% Arguments:    test case
% Returns:      none
%
% Description:  Tests that NumericSolver converges to steady
%               state for pure diffusion with Dirichlet BCs
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function TestGaussianQuadratureAccuracy(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function:     TestGaussianQuadratureAccuracy()
%
% Arguments:    test case
% Returns:      none
%
% Description:  Tests that Gaussian quadrature gives better
%               accuracy than trapezoidal for the same problem
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xmin = 0; xmax = 1;
    element_count = 5;
    order = 2; % quadratic elements
    
    mesh = Mesh(xmin, xmax, element_count, order, 1.0, 0.0);
    mesh.Generate();
    
    tmax = 0.1; dt = 0.01; theta = 0.5;
    
    lhs_bc = BoundaryCondition();
    lhs_bc.Type = BoundaryType.Dirichlet;
    lhs_bc.Value = 0.0;
    
    rhs_bc = BoundaryCondition();
    rhs_bc.Type = BoundaryType.Dirichlet;
    rhs_bc.Value = 1.0;
    
    % Solve with trapezoidal
    trap_method = IntegrationMethod();
    trap_method.type = IntegrationType.Trapezoidal;
    
    trap_solution = NumericSolver.SolveNumeric(mesh, tmax, dt, theta, ...
        lhs_bc, rhs_bc, @(x,t) 0, trap_method);
    
    % Solve with Gaussian
    gauss_method = IntegrationMethod();
    gauss_method.type = IntegrationType.Gaussian;
    gauss_method.gauss_points = 3;
    
    gauss_solution = NumericSolver.SolveNumeric(mesh, tmax, dt, theta, ...
        lhs_bc, rhs_bc, @(x,t) 0, gauss_method);
    
    % Get analytical solution
    analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);
    
    % Compute errors
    trap_error = L2Error(analytical_solution, trap_solution);
    gauss_error = L2Error(analytical_solution, gauss_solution);
    
    % Gaussian should be more accurate (lower final error)
    verifyLessThan(testCase, gauss_error.l2_error(end), ...
        trap_error.l2_error(end));
end

function TestQuadraticVsLinearElements(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function:     TestQuadraticVsLinearElements()
%
% Arguments:    test case
% Returns:      none
%
% Description:  Tests that quadratic elements give better
%               accuracy than linear elements for the same
%               problem
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xmin = 0; xmax = 1;
    element_count = 5;
    
    tmax = 0.1; dt = 0.01; theta = 0.5;
    
    lhs_bc = BoundaryCondition();
    lhs_bc.Type = BoundaryType.Dirichlet;
    lhs_bc.Value = 0.0;
    
    rhs_bc = BoundaryCondition();
    rhs_bc.Type = BoundaryType.Dirichlet;
    rhs_bc.Value = 1.0;
    
    integration_method = IntegrationMethod();
    integration_method.type = IntegrationType.Gaussian;
    integration_method.gauss_points = 3;
    
    % Linear elements
    mesh_linear = Mesh(xmin, xmax, element_count, 1, 1.0, 0.0);
    mesh_linear.Generate();
    
    solution_linear = NumericSolver.SolveNumeric(mesh_linear, tmax, dt, ...
        theta, lhs_bc, rhs_bc, @(x,t) 0, integration_method);
    
    analytical_linear = AnalyticalSolver.SolveAnalytical(mesh_linear, tmax, dt);
    error_linear = L2Error(analytical_linear, solution_linear);
    
    % Quadratic elements
    mesh_quadratic = Mesh(xmin, xmax, element_count, 2, 1.0, 0.0);
    mesh_quadratic.Generate();
    
    solution_quadratic = NumericSolver.SolveNumeric(mesh_quadratic, tmax, dt, ...
        theta, lhs_bc, rhs_bc, @(x,t) 0, integration_method);
    
    analytical_quadratic = AnalyticalSolver.SolveAnalytical(mesh_quadratic, tmax, dt);
    error_quadratic = L2Error(analytical_quadratic, solution_quadratic);
    
    % Quadratic should have lower error
    verifyLessThan(testCase, error_quadratic.l2_error(end), ...
        error_linear.l2_error(end));
end