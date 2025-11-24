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
    mesh = OneDimSimpleRefinedMeshGen(0, 1, 10);

    results = [];

    for t = 0:0.01:0.1
        % Solve transient diffusion equation using FEM
        timestep_results = [];

        for i = 1:mesh.ne
            timestep_results(i) = TransientAnalyticSoln((mesh.elem(i).x(1) + mesh.elem(i).x(2))/2, t);
        end

        results = [results; timestep_results];
    end

    % Display results as a heatmap
    figure;
    imagesc(results);
    colorbar;
    title('Transient Diffusion Solution Heatmap');
    xlabel('Element Index');
    ylabel('Time Step');

    saveas(gcf, 'main.fig');
    %saveas(gcf, 'cw1/report/resources/LaplaceEquationAnalyticalSolution.png');
    openfig('main.fig');

    dt = 0.01;
    steps = 10;
    theta = 0.5; % Crank-Nicholson
    D = 1;
    lambda = 0;

    solver = Solver(            ...
        dt,                     ...
        steps,                  ...
        theta,                  ...
        mesh,                   ...
        D,                      ...
        lambda                  ...
    );

    for step = 1:steps
        solver = solver.Update();
    end

    figure;
    imagesc(solver.solution');
    colorbar;
    title('Transient Diffusion FEM Solution Heatmap');
    xlabel('Time Step');
    ylabel('Global Node Index');
    saveas(gcf, 'fem_solution.fig');
    openfig('fem_solution.fig');
    

   fprintf('Exiting...\n');
end

function D = GetD(~)
    D = 1; % constant diffusion coefficient
end

function lambda = GetLambda(~)
    lambda = 0; % no reaction
end

function K = GetK(~)
    K = []; % placeholder
end

function M = GetM(~)
    M = []; % placeholder
end
