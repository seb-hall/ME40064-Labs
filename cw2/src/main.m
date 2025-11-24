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
    mesh = OneDimLinearMeshGen(0, 1, 50);



    results = [];

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

    t_vec = 0:dt:(steps*dt);
    x_vec = 1:mesh.ngn;

    % Display results as a heatmap
    figure;
    imagesc(t_vec, x_vec, results_analytic');
    colorbar;
    title('Transient Diffusion Solution Heatmap');
    xlabel('Time Step');
    ylabel('Element Index');
    caxis([0 1]); % Lock color scale

    saveas(gcf, 'main.fig');
    %saveas(gcf, 'cw1/report/resources/LaplaceEquationAnalyticalSolution.png');
    openfig('main.fig');

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
    imagesc(t_vec, x_vec, solver.solution);
    colorbar;
    title('Transient Diffusion FEM Solution Heatmap');
    xlabel('Time Step');
    ylabel('Element Index');
    caxis([0 1]); % Lock color scale
    saveas(gcf, 'fem_solution.fig');
    openfig('fem_solution.fig');

    % plot 
    

   fprintf('Exiting...\n');
end

function s = SourceFunction(x, t)
    s = 0; % No source term
end
