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
        function PlotHeatMap(solution, title_str, name)

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            % Plot a heat map of the solution values over time
            figure;
            imagesc(solution.time, solution.mesh.node_coords, solution.values);
            colorbar;
            xlabel("Time");
            ylabel("Position");
            caxis([0 1]) % lock color axis for consistency
            axis xy; % ensure y-axis is oriented correctly
            title(title_str);

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

        %% Plot solution at a specific position over time
        function PlotSampleOverTime(solution, x_sample, title_str, name)

            set(0, "DefaultAxesFontSize", 12);
            set(0, "DefaultTextFontSize", 12);

            % find x index
            x_index = (x_sample - solution.mesh.xmin) / (solution.mesh.xmax - solution.mesh.xmin) * solution.mesh.element_count + 1; % +1 for MATLAB indexing

            figure;
            plot_handle = plot(solution.time, solution.values(x_index, :));
            set(plot_handle, "LineWidth", 1.5);

            xlabel("Time (t)");

            
            ylabel("c(" + num2str(x_sample) + ", t)");
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