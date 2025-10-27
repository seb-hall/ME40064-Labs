function RunDiffusionReaction()
%% RunDiffusionReaction: solves a 1D mesh for the Diffusion/Reaction equation
%% as required for Q3.

    fprintf("Running Diffusion/Reaction Solver...\n");

    % array of different element counts
    elementCounts = [4, 8, 16];

    % iterate over each element count
    for Ne = elementCounts
        xmin = 0;
        xmax = 1;

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

        % clear current figure
        clf;

        % generate x vectors 
        x = linspace(xmin, xmax, Ne + 1);
        xreal = linspace(xmin, xmax, 500);
        
        % calculate real solution
        realSolution = exp(3)/(exp(6)-1)*(exp(3*xreal)-exp(-3*xreal));

        % plot solver and real solutions
        solverSolution = plot(x, solution);
        hold on;
        realHandle = plot(xreal, realSolution);

        % set chart title and axes
        title('Diffusion-Reaction Solution (' + string(Ne) + ' Elements)');
        xlabel('x');
        ylabel('c(x)');
        grid on;
        set(gcf, 'Position', [0, 0, 500, 350]);

        % set line properties
        set(realHandle, 'LineWidth', 1.5);
        set(realHandle, 'Color', [1, 0, 0, 0.5]);
        set(solverSolution, 'LineWidth', 1.5);
        set(solverSolution, 'Color', [0, 0, 1, 0.5]);

        % save and open figure
        saveas(gcf, 'DiffusionReactionSolution' + string(Ne) + 'Elements.fig');
        saveas(gcf, 'cw1/report/resources/DiffusionReactionSolution'...
            + string(Ne) + 'Elements.png');
        openfig('DiffusionReactionSolution' + string(Ne) + 'Elements.fig');
    end

end