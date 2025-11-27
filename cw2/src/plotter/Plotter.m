%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ME40064 Coursework 2
%
% File         :  Plotter.m
% Author       :  samh25
% Created      :  2025-11-26 (YYYY-MM-DD)
% License      :  MIT
% Description  :  A collection of static methods for plotting
%                 results for the coursework.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Plotter

    methods (Static)

        %% Plot entire solution as a heatmap - time as x-axis, position as y-axis and solution value as color
        function PlotHeatMap(solution, title_str, name, c_max)

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

        %% Plot full solution at specified time samples
        function PlotTimeSamples(solution, dt, time_samples, title_str, name)

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            figure;
            plot_handle = 0;

            for i = 1:length(time_samples)
                t_sample = time_samples(i);

                step_index = round(t_sample / dt) + 1; % +1 for MATLAB indexing

                plot_handle = plot(solution.mesh.node_coords, solution.values(:, step_index));
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

        %% Plot two solutions at a specific position over time
        function PlotSampleOverTime(solution_1, solution_2, x_sample, title_str, name, legend_strings)

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            % find x index
            x_index = round((x_sample - solution_1.mesh.xmin) / (solution_1.mesh.xmax - solution_1.mesh.xmin) * solution_1.mesh.element_count) + 1; % +1 for MATLAB indexing

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

    end
end