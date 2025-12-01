%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Plotter.m
% Author       :  11973
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A collection of static methods for plotting
%                 results for the coursework.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Plotter

    methods (Static)

        function PlotHeatMap(solution, title_str, name, c_max)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotHeatMap()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots a 2D heat map of solution values
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            % Plot a heat map of the solution values over time
            figure;
            imagesc(solution.time, solution.mesh.node_coords, solution.values);
            colorbar;
            xlabel("Time (t)");
            ylabel("Displacement (x)");
            caxis([0 c_max]) % lock color axis for consistency
            axis xy; % ensure y-axis is oriented correctly
            title(title_str);
            grid off;

            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");
            
        end

        function PlotTimeSamples(solution, dt, time_samples, title_str, name)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotTimeSamples()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots a solution at multiple time samples
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;
            plot_handle = 0;

            for i = 1:length(time_samples)
                t_sample = time_samples(i);

                step_index = round(t_sample / dt) + 1; % +1 for MATLAB indexing

                plot_handle = plot(solution.mesh.node_coords, ... 
                    solution.values(:, step_index));
                set(plot_handle, "LineWidth", 1.5);

                hold on;
            end

            xlabel("Position (x)");
            ylabel("c(x, t)");
            title(title_str);

            grid on;

            legend_strings = cell(1, length(time_samples));
            for i = 1:length(time_samples)
                legend_strings{i} = ['t = ', num2str(time_samples(i))];
            end

            legend(legend_strings, "Location", "northwest");

            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");

        end

        function PlotSampleOverTime(solution_1, solution_2, x_sample, title_str, name, legend_strings)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotSampleOverTime()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots two solutions at a given spatial sample over time
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            % find x index, +1 for MATLAB indexing
            x_index = round((x_sample - solution_1.mesh.xmin) / (solution_1.mesh.xmax ...
                - solution_1.mesh.xmin) * solution_1.mesh.element_count) + 1; 

            figure;
            plot_handle = plot(solution_1.time, solution_1.values(x_index, :));
            set(plot_handle, "LineWidth", 1.5);

            hold on;

            plot_handle = plot(solution_2.time, solution_2.values(x_index, :));
            set(plot_handle, "LineWidth", 1.5);

            xlabel("Time (t)");
            
            ylabel("c(" + num2str(x_sample) + ", t)");
            title(title_str);

            grid on;

            legend(legend_strings, "Location", "southeast");
            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");

        end

        function PlotConvergenceError(x_values, y_values, title_str, name, x_label, y_label)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotConvergenceError()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots convergence error on a log-log scale
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;

            plot_handle = loglog(x_values, y_values);
            set(plot_handle, "LineWidth", 1.5);
            
            xlabel(x_label);
            ylabel(y_label);
            title(title_str);
            grid on;

            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");

        end

        function PlotL2Errors(l2_errors, title_str, name, legend_strings)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotL2Errors()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots L2 errors over time for multiple simulations
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);  

            figure;
            
            for i = 1:length(l2_errors)
                l2_error = l2_errors(i);
                plot_handle = plot(l2_error.time, l2_error.l2_error);
                set(plot_handle, "LineWidth", 1.5);
                hold on;
            end

            xlabel("Time (t)");
            ylabel("L2 Error");
            title(title_str);

            grid on;

            legend(legend_strings, "Location", "northeast");
            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");

            openfig(name + ".fig");

        end

        function PlotTwoConvergenceLines(x_values, y1_values, y2_values, title_str, name, x_label, y_label, legend_strings)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotTwoConvergenceLines()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots two convergence lines on a log-log scale
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;
            
            loglog(x_values, y1_values, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
            hold on;
            loglog(x_values, y2_values, '-s', 'LineWidth', 1.5, 'MarkerSize', 8);
            
            xlabel(x_label);
            ylabel(y_label);
            title(title_str);
            legend(legend_strings, 'Location', 'best');
            grid on;

            set(gcf, 'Position', [0, 0, 500, 350]);

            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");
        end

        function PlotDoseEffectiveness(dose_values, kappa_values, title_str, name, x_label, y_label)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotDoseEffectiveness()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots dose vs effectiveness
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;

            plot_handle = plot(dose_values, kappa_values, "-o", 'LineWidth', 1.5);
            xlabel(x_label);
            ylabel(y_label);
            title(title_str);
            grid on;

            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");
        end

        function PlotSensitivityAnalysis(x_values, y_values, title_str, name, x_label, y_label, legend_strings)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotSensitivityAnalysis()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots sensitivity analysis results
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;

            % y_values is a matrix where each row is a different series
            for i = 1:size(y_values, 1)
                plot_handle = plot(x_values, y_values(i, :));
                set(plot_handle, "LineWidth", 1.5);
                hold on;
            end

            xlabel(x_label);
            ylabel(y_label);
            title(title_str);
            legend(legend_strings, 'Location', 'best');
            grid on;

            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");
        end
        
        function PlotKappaValues(x_values, y_values, title_str, name, x_label, y_label)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Function:     PlotKappaValues()
        %
        % Arguments:    Plotting parameters
        % Returns:      None
        %
        % Description:  Plots kappa values
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;
            
            plot_handle = plot(x_values, y_values, "-o");
            set(plot_handle, "LineWidth", 1.5);
            hold on;
    
            xlabel(x_label);
            ylabel(y_label);
            title(title_str);
            grid on;

            set(gcf, 'Position', [0, 0, 500, 350]);

            % Save figure
            saveas(gcf, name, "png");
            saveas(gcf, name, "fig");
            openfig(name + ".fig");
        end

    end
end