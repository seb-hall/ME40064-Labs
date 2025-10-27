function RunLaplace()

    fprintf("Running Laplace Solver...\n");

    xmin = 0;
    xmax = 1;
    Ne = 4;

    % create 1D uniform mesh
    mesh = OneDimLinearMeshGen(xmin, xmax, Ne);

    % set parameters for Laplace equation
    lambda = 0;
    D = 1; 

    leftBoundary = BoundaryCondition;
    leftBoundary.Type = BoundaryType.Neumann;
    leftBoundary.Value = 2; % dc/dx(c=0) = 2 -> flux at left boundary

    rightBoundary = BoundaryCondition;
    rightBoundary.Type = BoundaryType.Dirichlet;
    rightBoundary.Value = 0; % c(1) = 0 -> concentration at right boundary

    % call diffusion reaction solver
    solution = DiffusionReactionSolver(mesh, D, lambda, leftBoundary, rightBoundary);

    % plot solutions
    x = linspace(xmin, xmax, Ne + 1);
    plot(x, solution);
    title('Laplace Equation Solution');
    xlabel('x');
    ylabel('Solution c(x)');
    grid on;
    saveas(gcf, 'LaplaceEquationSolution.fig');
    openfig('LaplaceEquationSolution.fig');

end