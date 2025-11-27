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
            theta = 0.5;

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
            Plotter.PlotHeatMap(numeric_solution, "FEM Solution Heatmap", ...
                'cw2/report/resources/part1/NumericHeatmap', c_max);
            Plotter.PlotHeatMap(analytical_solution, "Analytical Solution Heatmap", ...
                'cw2/report/resources/part1/AnalyticalHeatmap', c_max);

            % plot solution samples at specified times
            sample_times = [0.05, 0.1, 0.3, 1.0];
            Plotter.PlotTimeSamples(numeric_solution, dt, sample_times, "FEM Solution Samples", ...
                'cw2/report/resources/part1/NumericSamples');
            Plotter.PlotTimeSamples(analytical_solution, dt, sample_times, "Analytical Solution Samples", ...
                'cw2/report/resources/part1/AnalyticalSamples');

            % plot both solutions at a specific position over time
            sample_x = 0.8;
            legend_strings = {'Numeric Solution', 'Analytical Solution'};
            Plotter.PlotSampleOverTime(numeric_solution, analytical_solution, ...
                sample_x, "Both Solutions at x = 0.8", 'cw2/report/resources/part1/BothX08', legend_strings);

        end

        function Part1Convergence()
            
            % time parameters
            tmax = 1.0;

            % mesh parameters
            xmin = 0.0;
            xmax = 1.0;
            element_count = 50;
            order = 1;
            
            % Crank-Nicholson method
            theta = 0.5;

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

            % calculate RMS error with varying mesh sizes and time steps

            element_counts = [5, 10, 20, 25, 50];
            time_steps = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005];

            num_cases = length(element_counts) * length(time_steps);
            rms_errr_table_elem_count = zeros(num_cases, 4);   % columns: elem_count, dt, dx, RMS error
            rms_errr_table_time_step = zeros(num_cases, 4);   % columns: elem_count, dt, dx, RMS error

            k = 1;

            for i = 1:length(element_counts)
                elem_count = element_counts(i);
                dt = 0.0005;

                % generate mesh
                mesh = Mesh(xmin, xmax, elem_count, order, D, lambda);
                mesh.Generate();

                % solve numerically
                numeric_solution = NumericSolver.SolveNumeric(...
                    mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @SourceFunction, integration_method);

                % solve analytically
                analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

                % compute RMS error
                [~, final_time] = min(abs(analytical_solution.time - tmax));

                c_numeric = numeric_solution.values(:, final_time);
                c_analytical = analytical_solution.values(:, final_time);

                error = c_numeric - c_analytical;
                rms_error = sqrt(mean(error.^2));

                rms_errr_table_elem_count(k,:) = [elem_count, dt, (xmax-xmin)/elem_count, rms_error];
                k = k + 1;

                fprintf('Elements: %d, dt: %.4f, dx: %.4f, RMS Error: %.6f\n', ...
                    elem_count, dt, (xmax-xmin)/elem_count, rms_error);
            end

            k = 1;

            if false
            for j = 1:length(time_steps)

                elem_count = 1 / 0.01;
                dt = time_steps(j);

                % generate mesh
                mesh = Mesh(xmin, xmax, elem_count, order, D, lambda);
                mesh.Generate();

                % solve numerically
                numeric_solution = NumericSolver.SolveNumeric(...
                    mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @SourceFunction, integration_method);

                % solve analytically
                analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

                % compute RMS error
                [~, final_time] = min(abs(analytical_solution.time - tmax));

                c_numeric = numeric_solution.values(:, final_time);
                c_analytical = analytical_solution.values(:, final_time);

                error = c_numeric - c_analytical;
                rms_error = sqrt(mean(error.^2));

                rms_errr_table_time_step(k,:) = [elem_count, dt, (xmax-xmin)/elem_count, rms_error];
                k = k + 1;

                fprintf('Elements: %d, dt: %.4f, dx: %.4f, RMS Error: %.6f\n', ...
                    elem_count, dt, (xmax-xmin)/elem_count, rms_error);
            end

            end 
            % plot element counts at min time step
            dx = rms_errr_table_elem_count(:, 3);  % element size
            err_spatial = rms_errr_table_elem_count(:, 4);

            Plotter.PlotConvergenceError(dx, err_spatial, ...
                "RMS Error vs Element Size (dt = 0.0005)", ...
                'cw2/report/resources/part1/ElementSizeConvergence', 'Element Size (x)', 'RMS Error at t = 1s');

            return;

            dt_vals = rms_errr_table_time_step(:, 2);  % time steps
            err_temporal = rms_errr_table_time_step(:, 4);

            Plotter.PlotConvergenceError(dt_vals, err_temporal, ...
                "RMS Error vs Time Step (dx = 0.01)", ...
                'cw2/report/resources/part1/TimeStepConvergence', 'Time Step (s)', 'RMS Error at t = 1s');
        end
    end

end


function s = SourceFunction(x, t)
    s = 0; % No source term
end
