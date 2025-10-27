function RunLaplace()
%% RunLaplace: solves a 1D mesh for Laplace's equation as required for Q2.

    fprintf("Running Laplace Solver...\n");

    % create 1D uniform mesh
    xmin = 0;
    xmax = 1;
    Ne = 4;
    mesh = OneDimLinearMeshGen(xmin, xmax, Ne);

    % set parameters for Laplace equation
    lambda = 0;
    D = 1; 

    % set boundary conditions
    leftBoundary = BoundaryCondition;
    leftBoundary.Type = BoundaryType.Neumann;
    leftBoundary.Value = 2; % dc/dx(c=0) = 2

    rightBoundary = BoundaryCondition;
    rightBoundary.Type = BoundaryType.Dirichlet;
    rightBoundary.Value = 0; % c(1) = 0

    % run diffusion reaction solver
    solution = DiffusionReactionSolver(mesh, D, lambda, leftBoundary, rightBoundary);

    % plot solution (assume uniform 1D mesh)
    x = linspace(xmin, xmax, Ne + 1);
    plotHandle = plot(x, solution);

    % set chart title and axes
    title('Laplace Equation Solution');
    xlabel('x');
    ylabel('c(x)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);
    set(plotHandle, 'LineWidth', 1.5);
    set(plotHandle, 'Color', [0, 0, 1]);

    % save and open figure
    saveas(gcf, 'LaplaceEquationSolution.fig');
    saveas(gcf, 'cw1/report/resources/LaplaceEquationSolution.png');
    openfig('LaplaceEquationSolution.fig');

    % plot real solution for comparison
    clf;
    realSolution = 2*(x - 1);
    plotHandle = plot(x, realSolution);

    % set chart title and axes
    title('Laplace Equation Analytical Solution');
    xlabel('x');
    ylabel('c(x)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);
    set(plotHandle, 'LineWidth', 1.5);
    set(plotHandle, 'Color', [1, 0, 0]);

    % save and open figure
    saveas(gcf, 'LaplaceEquationAnalyticalSolution.fig');
    saveas(gcf, 'cw1/report/resources/LaplaceEquationAnalyticalSolution.png');
    openfig('LaplaceEquationAnalyticalSolution.fig');

end