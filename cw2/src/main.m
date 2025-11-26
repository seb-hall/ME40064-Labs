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
    xmax = 1;
    element_count = 50;
    order = 1;

    theta = 0.5; % Crank-Nicholson
    D = 1;
    lambda = 0;

    mesh = Mesh(xmin, xmax, element_count, order, D, lambda);

    tmax = 1.0;
    dt = 0.01;

    analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

    Plotter.PlotHeatMap(analytical_solution, "Analytical Solution Heatmap", 'cw2/report/resources/AnalyticalHeatmap');

    lhs_boundary = BoundaryCondition();
    lhs_boundary.Type = BoundaryType.Dirichlet;
    lhs_boundary.Value = 0;

    rhs_boundary = BoundaryCondition();
    rhs_boundary.Type = BoundaryType.Dirichlet;
    rhs_boundary.Value = 1;

    numeric_solution = NumericSolver.SolveAnalytical(...
        mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @SourceFunction);

    Plotter.PlotHeatMap(numeric_solution, "Numeric Solution Heatmap", 'cw2/report/resources/NumericHeatmap');

    if true
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
