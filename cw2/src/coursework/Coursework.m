%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Coursework.m
% Author       :  samh25
% Created      :  2025-11-27 (YYYY-MM-DD)
% License      :  MIT
% Description  :  Static methods for each part of the coursework.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Coursework

    methods (Static)

        function Part1()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Part1()
        %
        % Arguments:    None
        % Returns:      None
        %
        % Description:  Runs Part 1 of the coursework, generating a simple
        %               mesh, running numeric and analytical solvers,
        %               and plotting the results.
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % time parameters
            tmax = 1.0;
            dt = 0.01;

            % mesh parameters
            xmin = 0.0;
            xmax = 1.0;
            element_count = 50;
            order = 1;
            
            % Crank-Nicholson method
            theta = 0.5 

            % diffusion and reaction coefficients
            D = 1.0; 
            lambda = 0.0;

            % concentrations
            c_max = 1.0;
            c_min = 0.0;

            % generate mesh
            mesh = Mesh(xmin, xmax, element_count, order, D, lambda);
            mesh.Generate();

            % solver parameters
            lhs_boundary = BoundaryCondition();
            lhs_boundary.Type = BoundaryType.Dirichlet;
            lhs_boundary.Value = c_min;

            rhs_boundary = BoundaryCondition();
            rhs_boundary.Type = BoundaryType.Dirichlet;
            rhs_boundary.Value = c_max;

            integration_method = IntegrationMethod();
            integration_method.type = IntegrationType.Trapezoidal;
            integration_method.gauss_points = 0; % not used for trapezoidal

            % solve numerically
            numeric_solution = NumericSolver.SolveNumeric(...
                mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @SourceFunction, integration_method);

            % solve analytically
            analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);
            
            % plot solutions as a heatmaps
            Plotter.PlotHeatMap(numeric_solution, "FEM Solution Heatmap", 'cw2/report/resources/part1/NumericHeatmap', c_max);
            Plotter.PlotHeatMap(analytical_solution, "Analytical Solution Heatmap", 'cw2/report/resources/part1/AnalyticalHeatmap', c_max);

            % sample times to plot
            sample_times = [0.05, 0.1, 0.3, 1.0];
            Plotter.PlotTimeSamples(numeric_solution, dt, sample_times, "FEM Solution Samples", 'cw2/report/resources/part1/NumericSamples');
            Plotter.PlotTimeSamples(analytical_solution, dt, sample_times, "Analytical Solution Samples", 'cw2/report/resources/part1/AnalyticalSamples');

            sample_x = 0.8;
            legend_strings = {'Numeric Solution', 'Analytical Solution'};
            Plotter.PlotSampleOverTime(numeric_solution, analytical_solution, sample_x, "Both Solutions at x = 0.8", 'cw2/report/resources/part1/BothX08', legend_strings);
            

        end
    end

end


function s = SourceFunction(x, t)
    s = 0; % No source term
end
