function RunLaplace()

    fprintf("Running Laplace Solver...\n");

    
    xmin = 0;
    xmax = 1;
    Ne = 4;

    % call uniform laplace solver
    solution = UniformLaplaceSolver(xmin, xmax, Ne);

    % plot solution
    x = linspace(xmin, xmax, Ne + 1);
    plot(x, solution);
    title('Laplace Equation Solution');
    xlabel('x');
    ylabel('Solution c(x)');
    grid on;
    saveas(gcf, 'LaplaceEquationSolution.fig');
    openfig('LaplaceEquationSolution.fig');

end