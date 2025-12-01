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

        function Part1Plots()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Part1Plots()
        %
        % Arguments:    None
        % Returns:      None
        %
        % Description:  Runs the start Part 1 of the coursework,  
        %               generating a simple mesh, running numeric and  
        %               analytical solvers, and plotting the results.
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
                mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, integration_method);

            % solve analytically
            analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);
            
            % plot solutions as a heatmaps
            Plotter.PlotHeatMap(numeric_solution, "FEM Solution Heatmap", ...
                "cw2/report/resources/part1/NumericHeatmap", c_max);
            Plotter.PlotHeatMap(analytical_solution, "Analytical Solution Heatmap", ...
                "cw2/report/resources/part1/AnalyticalHeatmap", c_max);

            % plot solution samples at specified times
            sample_times = [0.05, 0.1, 0.3, 1.0];
            Plotter.PlotTimeSamples(numeric_solution, dt, sample_times, "FEM Solution Samples", ...
                "cw2/report/resources/part1/NumericSamples");
            Plotter.PlotTimeSamples(analytical_solution, dt, sample_times, "Analytical Solution Samples", ...
                "cw2/report/resources/part1/AnalyticalSamples");

            % plot both solutions at a specific position over time
            sample_x = 0.8;
            legend_strings = {"Numeric Solution", "Analytical Solution"};
            Plotter.PlotSampleOverTime(numeric_solution, analytical_solution, ...
                sample_x, "Both Solutions at x = 0.8", "cw2/report/resources/part1/BothX08", legend_strings);

        end

        function Part1Convergence()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Part1Convergence()
        %
        % Arguments:    None
        % Returns:      None
        %
        % Description:  Runs a convergence study for Part 1 of the 
        %               coursework, calculating RMS error between numeric
        %               and analytical solutions over a range of element
        %               counts and time steps.
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            % vary element count with fixed time step
            for i = 1:length(element_counts)
                elem_count = element_counts(i);
                dt = 0.0005;

                % generate mesh
                mesh = Mesh(xmin, xmax, elem_count, order, D, lambda);
                mesh.Generate();

                % solve numerically
                numeric_solution = NumericSolver.SolveNumeric(...
                    mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, integration_method);

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

                fprintf("Elements: %d, dt: %.4f, dx: %.4f, RMS Error: %.6f\n", ...
                    elem_count, dt, (xmax-xmin)/elem_count, rms_error);
            end

            k = 1;

            % vary time step with fixed element count
            for j = 1:length(time_steps)

                elem_count = 1 / 0.01;
                dt = time_steps(j);

                % generate mesh
                mesh = Mesh(xmin, xmax, elem_count, order, D, lambda);
                mesh.Generate();

                % solve numerically
                numeric_solution = NumericSolver.SolveNumeric(...
                    mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, integration_method);

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

                fprintf("Elements: %d, dt: %.4f, dx: %.4f, RMS Error: %.6f\n", ...
                    elem_count, dt, (xmax-xmin)/elem_count, rms_error);

            end

            % plot element counts
            dx = rms_errr_table_elem_count(:, 3);  % element size
            err_spatial = rms_errr_table_elem_count(:, 4);

            Plotter.PlotConvergenceError(dx, err_spatial, ...
                "RMS Error vs Element Size (dt = 0.0005)", ...
                "cw2/report/resources/part1/ElementSizeConvergence", "Element Size (x)", "RMS Error at t = 1s");

            % plot time steps
            dt_vals = rms_errr_table_time_step(:, 2);  % time steps
            err_temporal = rms_errr_table_time_step(:, 4);

            Plotter.PlotConvergenceError(dt_vals, err_temporal, ...
                "RMS Error vs Time Step (dx = 0.01)", ...
                "cw2/report/resources/part1/TimeStepConvergence", "Time Step (s)", "RMS Error at t = 1s");
        end

        function Part2TimeIntegrationComparison()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Part2TimeIntegrationComparison()
        %
        % Arguments:    None
        % Returns:      None
        %
        % Description:  Runs a study comparing different time integration
        %               methods for Part 2 of the coursework.
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % time parameters
            tmax = 0.002;
            dt = 0.0001;

            % mesh parameters
            xmin = 0.0;
            xmax = 1.0;
            element_count = 10;
            order = 1;

            % diffusion and reaction coefficients
            D = 1.0; 
            lambda = 0.0;

            % concentrations
            c_max = 1.0;
            c_min = 0.0;

             % generate mesh
            mesh = Mesh(xmin, xmax, element_count, order, D, lambda);
            mesh.Generate();

            % solve analytically
            analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

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
            

            l2_errors = [];

            thetas = [0.0, 1.0, 0.5]; % Explicit Euler, Implicit Euler, Crank-Nicholson
            method_names = {"Forward Euler", "Crank-Nicolson", "Backward Euler"};

            for i = 1:length(thetas)
                theta = thetas(i);

                % solve numerically
                numeric_solution = NumericSolver.SolveNumeric(...
                    mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, integration_method);

                % compute L2 error
                l2_error = L2Error(analytical_solution, numeric_solution);
                l2_errors = [l2_errors, l2_error];
            end

            Plotter.PlotL2Errors(l2_errors, "L2 Error over Time", ...
                "cw2/report/resources/part2/L2ErrorTimeIntegration", ...
                method_names);

            % perform stability analysis 

            tmax = 1.0;
            element_count = 50;
            dt_list = [0.0001, 0.001, 0.01, 0.1, 0.25];

            % generate mesh
            mesh = Mesh(xmin, xmax, element_count, order, D, lambda);
            mesh.Generate();

            l2_errors_stability = [];

            for i = 1:length(thetas)
                theta = thetas(i);

                l2_errors_dt = [];

                for j = 1:length(dt_list)
                    dt = dt_list(j);

                    try

                        fprintf("Testing %s with dt = %.4f...\n", method_names{i}, dt);

                        numeric_solution = NumericSolver.SolveNumeric(...
                            mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, integration_method);

                            % CHECK FOR NaN/Inf at each timestep
                        for t_idx = 1:length(numeric_solution.time)
                            vals = numeric_solution.values(:, t_idx);
                            if any(isnan(vals))
                                fprintf("%s: NaN at step %d (t=%.4f)\n", method_names{i}, t_idx, numeric_solution.time(t_idx));
                                break;
                            end
                            if any(isinf(vals))
                                fprintf("%s: Inf at step %d (t=%.4f)\n", method_names{i}, t_idx, numeric_solution.time(t_idx));
                                break;
                            end
                            
                            if max(abs(vals)) > 1e10
                                fprintf("%s: Explosion at step %d (t=%.4f), max=%.2e\n", method_names{i}, t_idx, numeric_solution.time(t_idx), max(abs(vals)));
                                break;
                            end
                        end

                        analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

                        l2_error = L2Error(analytical_solution, numeric_solution);
                        
                        l2_errors_dt = [l2_errors_dt, l2_error];
                    catch
                        l2_errors_dt = [l2_errors_dt, NaN];
                        fprintf("%s EXPLODED \n", method_names{i});
                    end
                end

                l2_errors_stability = [l2_errors_stability; l2_errors_dt];
            end

        end

        function Part2GaussianQuadrature()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Part2GaussianQuadrature()
        %
        % Arguments:    None
        % Returns:      None
        %
        % Description:  Runs a study comparing L2 error with and without
        %               Gaussian Quadrature.
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % time parameters
            tmax = 0.001;
            dt = 0.0001;

            % mesh parameters
            xmin = 0.0;
            xmax = 1.0;
            element_count = 5;
            order = 2;

            % diffusion and reaction coefficients
            D = 0.5; 
            lambda = 0.0;

            % concentrations
            c_max = 1.0;
            c_min = 0.0;

            % generate mesh
            mesh = Mesh(xmin, xmax, element_count, order, D, lambda);
            mesh.Generate();

            % solve analytically
            analytical_solution = AnalyticalSolver.SolveAnalytical(mesh, tmax, dt);

            % solver parameters

            theta = 0.5; % Crank-Nicholson

            lhs_boundary = BoundaryCondition();
            lhs_boundary.Type = BoundaryType.Dirichlet;
            lhs_boundary.Value = c_min;

            rhs_boundary = BoundaryCondition();
            rhs_boundary.Type = BoundaryType.Dirichlet;
            rhs_boundary.Value = c_max;

            % trapezoidal method
            trapezoidal_method = IntegrationMethod();
            trapezoidal_method.type = IntegrationType.Trapezoidal;
            trapezoidal_method.gauss_points = 0; % not used for trapezoidal

            trapezoidal_solution = NumericSolver.SolveNumeric(...
                mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, trapezoidal_method);

            % gaussian quadrature method
            gaussian_method = IntegrationMethod();
            gaussian_method.type = IntegrationType.Gaussian;
            gaussian_method.gauss_points = 3; % 3-point Gaussian quadrature

            gaussian_solution = NumericSolver.SolveNumeric(...
                mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, gaussian_method);

            % compute L2 error
            l2_error_trapezoidal = L2Error(analytical_solution, trapezoidal_solution);
            l2_error_gaussian = L2Error(analytical_solution, gaussian_solution);
            l2_errors = [l2_error_trapezoidal, l2_error_gaussian];

            method_names = {"2-point Trapezoidal Integration", "3-point Gaussian Quadrature"};

            Plotter.PlotL2Errors(l2_errors, "Gaussian vs Trapezoidal Integration", ...
                "cw2/report/resources/part2/L2ErrorGaussianTrapezoidal", ...
                method_names);
        end

        function Part3InitialResults()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     Part2GaussianQuadrature()
        %
        % Arguments:    None
        % Returns:      None
        %
        % Description:  Runs a study comparing L2 error with and without
        %               Gaussian Quadrature.
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % Generate mesh
            xmin = 0;
            xmax = 0.01;
            element_count = 50;
            order = 2;

            theta = 0.5; % Crank-Nicholson
            D = 1;
            lambda = 0;

            epidermis_layer = MeshLayer(0.0, 4e-6, 0.0, 0.02, 1.0);
            dermis_layer = MeshLayer(0.00166667, 5e-6, 0.01, 0.02, 1.0);
            sub_cutaneous_layer = MeshLayer(0.005, 2e-6, 0.01, 0.02, 1.0);

            layers = [epidermis_layer, dermis_layer, sub_cutaneous_layer];

            mesh = MultilayerMesh(xmin, xmax, element_count, order, D, lambda, layers);
            mesh.Generate();

            tmax = 30.0;
            dt = 0.01; % works well with element_count = 50

            % concentrations
            c_max = 30.0;
            c_min = 0.0;

            lhs_boundary = BoundaryCondition();
            lhs_boundary.Type = BoundaryType.Dirichlet;
            lhs_boundary.Value = c_max;

            rhs_boundary = BoundaryCondition();
            rhs_boundary.Type = BoundaryType.Dirichlet;
            rhs_boundary.Value = c_min;

            integration_method = IntegrationMethod();
            integration_method.type = IntegrationType.Gaussian;
            integration_method.gauss_points = order + 1;

            numeric_solution = NumericSolver.SolveNumeric(...
                mesh, tmax, dt, theta, lhs_boundary, rhs_boundary, @(~, ~) 0.0, integration_method);


            kappa = DoseEvaluator.EvaluateSolution(numeric_solution, 0.005, 4.0, dt);
                
            fprintf('Kappa: %.2f\n', kappa);

            Plotter.PlotHeatMap(numeric_solution, "Drug Concentration Heatmap", 'cw2/report/resources/part3/InitialNumericHeatmap', c_max);
        end
            

    end

end
