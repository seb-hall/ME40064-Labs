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
    [
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
    ],
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

== Results

== Source Code

= Q3: Verifying the FEM Solver for the Diffusion-Reaction Equation

== Results

== Source Code

// MARK: CONCLUSION
== Possible Improvements



// MARK: REFERENCES
= References

#bibliography(
    "resources/ME40064-References.yml",
    title: none,
    style: "ieee"
)