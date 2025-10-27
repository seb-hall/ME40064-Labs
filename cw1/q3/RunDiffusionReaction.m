function RunDiffusionReaction()
%% RunDiffusionReaction: solves a 1D mesh for the Diffusion/Reaction equation equation as required for Q3.

    fprintf("Running Diffusion/Reaction Solver...\n");

    xmin = 0;
    xmax = 1;
    Ne = 40;

    % create 1D uniform mesh
    mesh = OneDimLinearMeshGen(xmin, xmax, Ne);

    % set parameters for Laplace equation
    lambda = -9;
    D = 1; 

    leftBoundary = BoundaryCondition;
    leftBoundary.Type = BoundaryType.Dirichlet;
    leftBoundary.Value = 0; % c(0) = 0 -> concentration at right boundary

    rightBoundary = BoundaryCondition;
    rightBoundary.Type = BoundaryType.Dirichlet;
    rightBoundary.Value = 1; % c(1) = 1 -> concentration at right boundary

    % call diffusion reaction solver
    solution = DiffusionReactionSolver(mesh, D, lambda, leftBoundary, rightBoundary);

    for i = 1:length(solution)
        fprintf("Node %d: x = %.4f, c = %.6f\n", i, mesh.nvec(i), solution(i));
    end

    % plot solution (assume uniform 1D mesh)
    x = linspace(xmin, xmax, Ne + 1);
    plot(x, solution);

    % set chart title and axes
    title('Diffusion/Reaction Solution');
    xlabel('x');
    ylabel('c(x)');
    grid on;
    set(gcf, 'Position', [0, 0, 500, 350]);

    % save and open figure
    saveas(gcf, 'DiffusionReactionSolution.fig');
    saveas(gcf, 'cw1/report/resources/DiffusionReactionSolution.png');
    openfig('DiffusionReactionSolution.fig');

end