#import "@preview/plotst:0.2.0": *

#import "resources/ME40064-References.yml"

#set page(
    paper: "a4",
    margin: (x: 1.25cm, top: 1.5cm, bottom: 1.5cm),
    columns: 1,
    header:  context {
        if(counter(page).get().at(0) != 1) [
            *ME40064 System Modelling and Simulation - Coursework 1*
            #h(1fr)
            #counter(page).display(
                "1/1",
                both: true,
            )
        ]
    }
)

#set text(
    size: 11pt
)

#set par(
    justify: true,
    leading: 0.52em,
)

#set heading(
    numbering: "1."
)

#place(
    top + center,
    float: true,
    scope: "parent",
    {
        text(17pt)[
            *ME40064 System Modelling and Simulation - Coursework 1* \
        ]

        text(13pt)[
            Seb Hall #link("mailto:samh25@bath.ac.uk"), 4th November 2025\
            Department of Mechanical Engineering, University of Bath \
        ]
    }
)

#show raw.where(block: true): block.with(fill: luma(240), inset: 1em, radius: 0.5em, width: 100%)
#show figure: set block(breakable: true)
#show link: underline

// MARK: INTRO
= Introduction

// MARK: EX 1
= Q1: Local Element Matrix Functions

== Q1a: MATLAB 2x2 Element Matrix for the Diffusion Operator


=== Unit Test Result
#figure(
    image("resources/q1a_test_result.png", width: 100%),
    caption: [Diffusion Operator Unit Test Results],
)  <1a-test-results>

=== Function Code

#figure(
    align(left, [
        ```matlab
        function matrix = DiffusionElemMatrix(D, eID, msh)
        %% DiffusionElemMatrix - calculates a 2x2 element matrix for the diffusion 
        %% operator in a 1d finite element mesh.
            
            % create base matrix
            matrix = [1, -1; -1, 1];    

            % calulate element size
            elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

            % apply matrix scaling
            matrix = matrix * (D / elemSize);
        
        end
        ```
    ]),
    caption: [Diffusion Operator Element Matrix Function],
    supplement: [Code Snippet]
)  <1a-code>

== Q1b: MATLAB 2x2 Element Matrix for the Linear Reaction Operator

=== Unit Test Result
#figure(
    image("resources/q1b_test_result.png", width: 100%),
    caption: [Linear Reaction Operator Unit Test Results], 
)  <1b-test-results>


=== Function Code

#figure(
    align(left, [
        ```matlab
        function matrix = ReactionElemMatrix(lambda, eID, msh)
        %% ReactionElemMatrix - calculates a 2x2 element matrix for the linear
        %% reaction operator in a 1d finite element mesh.
            
            % create base matrix
            matrix = [2, 1; 1, 2];    

            % calulate element size
            elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

            % apply matrix scaling
            matrix = matrix * (lambda * elemSize / 6);
        
        end
        ```
    ]),
    caption: [Linear Reaction Operator Element Matrix Function],
    supplement: [Code Snippet]
)  <1b-code>

=== Unit Test 1

#figure(
    align(left, [
        ```matlab
        %% Test 1: test symmetry of the matrix
        % % Test that this matrix is symmetric

        tol = 1e-14; % test tolerance 
        lambda = 2; % reaction coefficient
        eID = 1; % element ID

        xmin = 0;
        xmax = 1;
        Ne = 10;
        msh = OneDimLinearMeshGen(xmin, xmax, Ne);

        elemat = ReactionElemMatrix(lambda, eID, msh);

        assert(abs(elemat(1,2) - elemat(2,1)) <= tol)
        ```
    ]),
    caption: [Linear Reaction Operator Unit Test 1],
    supplement: [Code Snippet]
)  <1b-test-1>

=== Unit Test 2

#figure(
    align(left, [
        ```matlab
        %% Test 2: test 2 different elements of the same size produce same matrix
        % % Test that for two elements of an equispaced mesh, the element matrices 
        % % are calculated are the same.

        tol = 1e-14; % test tolerance
        lambda = 5; % reaction coefficient
        eID = 1; % element ID

        xmin = 0;
        xmax = 1;
        Ne = 10;
        msh = OneDimLinearMeshGen(xmin, xmax, Ne);

        elemat1 = ReactionElemMatrix(lambda, eID, msh);

        eID = 2; %element ID
        elemat2 = ReactionElemMatrix(lambda, eID, msh);

        diff = elemat1 - elemat2;
        diffnorm = sum(sum(diff.*diff));
        assert(abs(diffnorm) <= tol)
        ```
    ]),
    caption: [Linear Reaction Operator Unit Test 2],
    supplement: [Code Snippet]
)  <1b-test-2>


=== Unit Test 3

#figure(
    align(left, [
        ```matlab
        %% Test 3: test that one matrix is evaluated correctly
        % % Test that element 1 of the (equispaced) three element mesh 
        % % problem is evaluated correctly

        tol = 1e-14; % test tolerance 
        lambda = 2.5; % reaction coefficient
        eID = 1; % element ID

        xmin = 0;
        xmax = 1;
        Ne = 3;
        msh = OneDimLinearMeshGen(xmin, xmax, Ne);

        elemat1 = ReactionElemMatrix(lambda, eID, msh);

        matrix = [2, 1; 1, 2];
        elemSize = (xmax - xmin) / Ne;
        elemat2 = matrix * (lambda * elemSize / 6);

        diff = elemat1 - elemat2; % calculate the difference between the two matrices
        diffnorm = sum(sum(diff.*diff)); % calculate the total squared error between the matrices
        assert(abs(diffnorm) <= tol)
        ```
    ]),
    caption: [Linear Reaction Operator Unit Test 3],
    supplement: [Code Snippet]
)  <1b-test-3>

=== Unit Test 4

#figure(
    align(left, [
        ```matlab
        %% Test 4: test that different sized elements in a mesh are evaluated correctly - element 1
        % % Test that elements in a non-equally spaced mesh are evaluated correctly

        tol = 1e-14; % test tolerance
        lambda = 1; % reaction coefficient
        eID = 1; % element ID

        xmin = 0;
        xmax = 1;
        Ne = 5;
        msh = OneDimSimpleRefinedMeshGen(xmin, xmax, Ne);

        elemat1 = ReactionElemMatrix(lambda, eID, msh);

        elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1); % get element size from mesh
        elemat2 = [2, 1; 1, 2] * (lambda * elemSize / 6); % calculate resultant matrix

        diff = elemat1 - elemat2; % calculate the difference between the two matrices
        diffnorm = sum(sum(diff.*diff)); % calculate the total squared error between the matrices
        assert(abs(diffnorm) <= tol)
        ```
    ]),
    caption: [Linear Reaction Operator Unit Test 4],
    supplement: [Code Snippet]
)  <1b-test-4>

=== Unit Test 5

#figure(
    align(left, [
        ```matlab
        %% Test 5: test that different sized elements in a mesh are evaluated correctly - element 4
        % % Test that elements in a non-equally spaced mesh are evaluated correctly

        tol = 1e-14; % test tolerance
        lambda = 1; % reaction coefficient
        eID = 4; % element ID

        xmin = 0;
        xmax = 1;
        Ne = 5;
        msh = OneDimSimpleRefinedMeshGen(xmin, xmax, Ne);

        elemat1 = ReactionElemMatrix(lambda, eID, msh);

        elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1); % get element size from mesh
        elemat2 = [2, 1; 1, 2] * (lambda * elemSize / 6); % calculate resultant matrix

        diff = elemat1 - elemat2; % calculate the difference between the two matrices
        diffnorm = sum(sum(diff.*diff)); % calculate the total squared error between the matrices
        assert(abs(diffnorm) <= tol)
        ```
    ]),
    caption: [Linear Reaction Operator Unit Test 5],
    supplement: [Code Snippet]
)  <1b-test-5>

= Q2: Solving Laplace's Equation using FEM

This task was to solve Laplace's equation:

$
    (delta ^2c)/(delta x^2) = 0
$

With the following boundary conditions: 

$
    (delta c)/(delta x)(x = 0) = 2
$

$
    c(x = 1) = 0
$

Over a 1D uniform mesh with 4 elements.

== Results

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    align: horizon,
    [
        #figure(
            image("resources/LaplaceEquationSolution.png", width: 100%),
            caption: [Solution of the Laplace Equation \ using FEM Solver],

        )  <and-state-space>
    ],
    [
        #figure(
            image("resources/LaplaceEquationAnalyticalSolution.png", width: 100%),
            caption: [Solution of the Laplace Equation using \ an Analytical Method, *y = 2(x + 1)*],

        )  <and-state-space>
    ] 
    )
])

These results show that for this example, the solver is able to produce a perfect solution for the Laplace equation. This is due to the function being linear, and therefore perfectly representable by discrete sample points.

== Source Code

=== Solver

#figure(
    align(left, [
        ```matlab
            classdef BoundaryType
                enumeration
                        Dirichlet, Neumann
                end
            end
        ```
    ]),
    caption: [BoundaryType.m],
    supplement: [Code Snippet]
)  <boundary-type-m>

#figure(
    align(left, [
        ```matlab
            classdef BoundaryCondition
                properties
                    Type % Boundary condition type (Dirichlet or Neumann)
                    Value % Boundary condition value
                end
            end
        ```
    ]),
    caption: [BoundaryCondition.m],
    supplement: [Code Snippet]
)  <boundary-condition-m>

#figure(
    align(left, [
        ```matlab
            function solution = DiffusionReactionSolver(mesh, D, lambda, leftBoundary, rightBoundary)
            %% DiffusionReactionSolver - solves the steady-state diffusion-reaction equation for a 1D mesh
            % % Inputs:
            %   mesh - 1D finite element mesh structure
            %   D - diffusion coefficient
            %   lambda - reaction rate
            %   leftBoundary - left boundary condition (BoundaryCondition object)
            %   rightBoundary - right boundary condition (BoundaryCondition object)
            %
            % % Outputs:
            %   solution - solution vector at mesh nodes

                % calculate number of nodes and elements
                Ne = mesh.ne;
                Nn = Ne + 1;

                % initialise global matrix
                globalMatrix = zeros(Nn, Nn);

                % assemble global matrix
                for eID = 1:Ne

                    % get element matrices
                    diffusionElementMatrix = DiffusionElemMatrix(D, eID, mesh);
                    reactionElementMatrix = ReactionElemMatrix(lambda, eID, mesh);

                    % combine element matrices
                    elemMatrix = diffusionElementMatrix - reactionElementMatrix;
                    
                    % insert into global matrix
                    globalMatrix(eID, eID) = globalMatrix(eID, eID) + elemMatrix(1, 1);
                    globalMatrix(eID, eID + 1) = globalMatrix(eID, eID + 1) + elemMatrix(1, 2);
                    globalMatrix(eID + 1, eID) = globalMatrix(eID + 1, eID) + elemMatrix(2, 1);
                    globalMatrix(eID + 1, eID + 1) = globalMatrix(eID + 1, eID + 1) + elemMatrix(2, 2);

                end

                % initialise source vector
                sourceVector = zeros(Nn, 1);

                % here we would assemble the source vector if there were any source terms
                % however, for now we assume there are none, so it remains zero

                % Apply left boundary condition
                switch leftBoundary.Type

                    case BoundaryType.Neumann

                        % directly modify source vector for Neumann
                        sourceVector(1) = sourceVector(1) - leftBoundary.Value;

                    case BoundaryType.Dirichlet

                        % apply Dirichlet condition to source vector
                        for j = 2:Nn
                            sourceVector(j) = sourceVector(j) - globalMatrix(j, 1) * leftBoundary.Value;
                        end

                        % modify global matrix
                        globalMatrix(1, :) = 0;
                        globalMatrix(:, 1) = 0;
                        globalMatrix(1, 1) = 1;

                        % set value at first node
                        sourceVector(1) = leftBoundary.Value; 
                end

                % Apply right boundary condition
                switch rightBoundary.Type

                    case BoundaryType.Neumann

                        % directly modify source vector for Neumann
                        sourceVector(Nn) = sourceVector(Nn) - rightBoundary.Value;

                    case BoundaryType.Dirichlet

                        % apply Dirichlet condition to source vector
                        for j = 2:(Nn-1)
                            sourceVector(j) = sourceVector(j) - globalMatrix(j, Nn) * rightBoundary.Value;
                        end

                        % modify global matrix
                        globalMatrix(Nn, :) = 0;
                        globalMatrix(:, Nn) = 0;
                        globalMatrix(Nn, Nn) = 1;

                        % set value at last node
                        sourceVector(Nn) = rightBoundary.Value; 
                end

                % solve system of equations
                solution = globalMatrix \ sourceVector;
            end


        ```
    ]),
    caption: [DiffusionReactionSolver.m],
    supplement: [Code Snippet]
)  <diffusion-reaction-solver-m>

=== Implementation

#figure(
    align(left, [
        ```matlab
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
        ```
    ]),
    caption: [RunLaplace.m],
    supplement: [Code Snippet]
)  <run-laplace-m>

= Q3: Verifying the FEM Solver for the Diffusion-Reaction Equation

== Results



== Source Code

== Possible Improvements

= Use of Generative AI

This coursework was completed in Visual Studio Code (with the #link("https://marketplace.visualstudio.com/items?itemName=MathWorks.language-matlab", "MATLAB Extension")), using Typst for report writing. The #link("https://github.com/features/copilot", "GitHub Copilot") AI tool was enabled for this, and provided generative suggestions for code snippets and report phrasing.

// MARK: REFERENCES
= References

#bibliography(
    "resources/ME40064-References.yml",
    title: none,
    style: "ieee"
)