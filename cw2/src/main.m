%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  main.m
% Author       :  samh25
% Created      :  2025-11-24 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Main function for solving transient diffusion equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
   fprintf('ME40064 Solver...\n');

    % Generate mesh

    xmin = 0;
    xmax = 0.01;
    element_count = 50;
    order = 2;

    theta = 0.5; % Crank-Nicholson
    D = 1;
    lambda = 0;

    epidermis_layer = LayerProperties(0.0, 4e-6, 0.0, 0.02);
    dermis_layer = LayerProperties(0.00166667, 5e-6, 0.01, 0.02);
    sub_cutaneous_layer = LayerProperties(0.005, 2e-6, 0.01, 0.02);

    layer_properties = [epidermis_layer, dermis_layer, sub_cutaneous_layer];

    mesh = MultilayerMesh(xmin, xmax, element_count, order, D, lambda, layer_properties);
    mesh.Generate();

    tmax = 30.0;
    dt = 0.01; % works well with element_count = 50

    analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

    lhs_boundary = BoundaryCondition();
    lhs_boundary.Type = BoundaryType.Dirichlet;
    lhs_boundary.Value = 30;

    rhs_boundary = BoundaryCondition();
    rhs_boundary.Type = BoundaryType.Dirichlet;
    rhs_boundary.Value = 0;

    integration_method = IntegrationMethod();
    integration_method.type = IntegrationType.Gaussian;
    integration_method.gauss_points = order + 1;

    numeric_solution = NumericSolver.SolveNumeric(...
        mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @SourceFunction, integration_method);


    l2_error = L2Error(analytical_solution, numeric_solution);

    kappa = DoseEvaluator.EvaluateSolution(numeric_solution, 0.005, 4.0, dt);
        
    fprintf('Kappa: %.2f\n', kappa);


    %Plotter.PlotL2Error(l2_error, "L2 Error over Time", 'cw2/report/resources/L2Error');

    if true

        %Plotter.PlotHeatMap(analytical_solution, "Analytical Solution Heatmap", 'cw2/report/resources/AnalyticalHeatmap');
        Plotter.PlotHeatMap(numeric_solution, "Numeric Solution Heatmap", 'cw2/report/resources/NumericHeatmap');

    end

    if false

         % sample times to plot
        sample_times = [0.05, 0.1, 0.3, 1.0];
        Plotter.PlotTimeSamples(analytical_solution, dt, sample_times, "Analytical Solution Samples", 'cw2/report/resources/AnalyticalSamples');
        Plotter.PlotTimeSamples(numeric_solution, dt, sample_times, "Numeric Solution Samples", 'cw2/report/resources/NumericSamples');

        sample_x = 0.8;
        Plotter.PlotSampleOverTime(analytical_solution, sample_x, "Analytical Solution at x = 0.8", 'cw2/report/resources/AnalyticalX08');
        Plotter.PlotSampleOverTime(numeric_solution, sample_x, "Numeric Solution at x = 0.8", 'cw2/report/resources/NumericX08');
    end

    fprintf('Exiting...\n');
end

function s = SourceFunction(x, t)
    s = 0; % No source term
end
