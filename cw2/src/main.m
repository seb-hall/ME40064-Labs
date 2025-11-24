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

    x_size = (xmax - xmin) / element_count;
    x_vals = xmin:x_size:xmax;

    mesh = OneDimMeshGen(xmin, xmax, element_count, 2);

    dt = 0.01;
    tmax = 1.0;
    steps = tmax / dt;

    % (Assuming a simple mesh where elements are sequential)
    node_coords = zeros(1, mesh.ngn);
    node_coords(1) = mesh.elem(1).x(1);
    for k = 1:mesh.ne
        node_coords(k+1) = mesh.elem(k).x(2);
    end

    results_analytic = [];
    
    % Define time vector matching the solver settings
    time_vector = 0:dt:(steps*dt); 

    for t = time_vector
        timestep_results = zeros(1, mesh.ngn);
        
        % Loop over NODES, not elements
        for i = 1:mesh.ngn
            timestep_results(i) = TransientAnalyticSoln(node_coords(i), t);
        end

        results_analytic = [results_analytic; timestep_results];
    end

    set(0, "DefaultAxesFontSize", 12);
    set(0, "DefaultTextFontSize", 12);

    t_vec = 0:dt:(steps*dt);

    % Display results as a heatmap
    figure;
    imagesc(t_vec, x_vals, results_analytic');
    colorbar;
    title('Transient Diffusion Solution Heatmap');
    xlabel('Time (s)');
    ylabel('Position (x)');
    caxis([0 1]); % Lock color scale
    axis xy;
    set(gcf, 'Position', [0, 0, 500, 350]);

    saveas(gcf, 'cw2/report/resources/AnalyticalHeatmap.fig');
    saveas(gcf, 'cw2/report/resources/AnalyticalHeatmap.png');
    openfig('cw2/report/resources/AnalyticalHeatmap.fig');

    theta = 0.5; % Crank-Nicholson
    D = 1;
    lambda = 0;

    lhs_boundary = BoundaryCondition();
    lhs_boundary.Type = BoundaryType.Dirichlet;
    lhs_boundary.Value = 0;

    rhs_boundary = BoundaryCondition();
    rhs_boundary.Type = BoundaryType.Dirichlet;
    rhs_boundary.Value = 1;

    solver = Solver(            ...
        dt,                     ...
        uint64(steps),          ...
        theta,                  ...
        mesh,                   ...
        D,                      ...
        lambda,                 ...
        lhs_boundary,           ...
        rhs_boundary,           ...
        @SourceFunction         ...
    );

    for step = 1:steps
        solver = solver.Update();
    end

    figure;
    imagesc(t_vec, x_vals, solver.solution);
    colorbar;
    title('Transient Diffusion FEM Solution Heatmap');
    xlabel('Time Step');
    ylabel('Position (x)');
    caxis([0 1]); % Lock color scale
    axis xy;
    set(gcf, 'Position', [0, 0, 500, 350]);
    saveas(gcf, 'cw2/report/resources/SolverHeatmap.fig');
    saveas(gcf, 'cw2/report/resources/SolverHeatmap.png');
    openfig('cw2/report/resources/SolverHeatmap.fig');

    % sample times to plot
    sample_times = [0.05, 0.1, 0.3, 1.0];

    

    % plot analytical solutions at sample times
    figure;
    plotHandle = 0; 

    for i = 1:length(sample_times)
        t_sample = sample_times(i);
        step_index = round(t_sample / dt) + 1; % +1 for MATLAB indexing
        
        plotHandle = plot(x_vals, results_analytic(step_index, :));
        set(plotHandle, 'LineWidth', 1.5);

        hold on;
    end

     % set chart title and axes
    title("Analytical Solutions");
    xlabel('x');
    ylabel('c(x, t)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);

    legendStrings = cell(1, length(sample_times));
    for i = 1:length(sample_times)
        legendStrings{i} = ['t = ', num2str(sample_times(i))];
    end 

    legend(legendStrings, "Location", "northwest");

    % save and open figure
    saveas(gcf, 'cw2/report/resources/AnalyticalSamples.fig');
    saveas(gcf, 'cw2/report/resources/AnalyticalSamples.png');
    openfig('cw2/report/resources/AnalyticalSamples.fig');

    figure;
    plotHandle = 0;

    for i = 1:length(sample_times)
        step_index = round(sample_times(i) / dt) + 1; % +1 for MATLAB indexing
       
        plotHandle = plot(x_vals, solver.solution(:, step_index));
        set(plotHandle, 'LineWidth', 1.5);

        hold on;
    end

     % set chart title and axes
    title("FEM Solutions");
    xlabel('x');
    ylabel('c(x, t)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);

    legendStrings = cell(1, length(sample_times));
    for i = 1:length(sample_times)
        legendStrings{i} = ['t = ', num2str(sample_times(i))];
    end 

    legend(legendStrings, "Location", "northwest");

    % save and open figure
    saveas(gcf, 'cw2/report/resources/SolverSamples.fig');
    saveas(gcf, 'cw2/report/resources/SolverSamples.png');
    openfig('cw2/report/resources/SolverSamples.fig');

    % sample x = 0.8 over time

    x_index = (0.8 - xmin)/(xmax - xmin) * element_count + 1; % +1 for MATLAB indexing

    figure;
    plotHandle = plot(t_vec, results_analytic(:, x_index));
    set(plotHandle, 'LineWidth', 1.5);
    title('Analytical Solution at x = 0.8');
    xlabel('Time (s)');
    ylabel('c(0.8, t)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);
    saveas(gcf, 'cw2/report/resources/AnalyticalX08.fig');
    saveas(gcf, 'cw2/report/resources/AnalyticalX08.png');
    openfig('cw2/report/resources/AnalyticalX08.fig');

    figure;
    plotHandle = plot(t_vec, solver.solution(x_index, :));
    set(plotHandle, 'LineWidth', 1.5);
    title('FEM Solution at x = 0.8');
    xlabel('Time (s)');
    ylabel('c(0.8, t)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);
    saveas(gcf, 'cw2/report/resources/SolverX08.fig');
    saveas(gcf, 'cw2/report/resources/SolverX08.png');
    openfig('cw2/report/resources/SolverX08.fig');
    


    fprintf('Exiting...\n');
end

function s = SourceFunction(x, t)
    s = 0; % No source term
end
