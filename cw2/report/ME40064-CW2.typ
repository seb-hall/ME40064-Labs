#import "@preview/plotst:0.2.0": *

#import "@preview/wordometer:0.1.4" as wordometer
#import wordometer: word-count, total-words, total-characters

#show: word-count.with(exclude: (bibliography, figure, raw, table, outline, <executive-summary>))

#import "resources/ME40064-References.yml"

#set page(
    paper: "a4",
    margin: (x: 1.25cm, top: 1.5cm, bottom: 1.5cm),
    columns: 2,
    header:  context {
        if(counter(page).get().at(0) != 1) [
            *ME40064 System Modelling and Simulation - Coursework 2*
            #h(1fr)
            #counter(page).display(
                "1/1",
                both: true,
            )
        ]
    }
)

#let style-number(number) = text(gray)[#number]

#show raw.where(block: true): it => grid(
  columns: 2,
  align: (right, left),
  gutter: 0.5em,
  ..it.lines
    .enumerate()
    .map(((i, line)) => (style-number(i + 1), line))
    .flatten()
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
            *ME40064 System Modelling and Simulation - Coursework 2* \
        ]

        text(13pt)[
            #total-words Words, Candidate No. 11973, 2nd December 2025\
            Department of Mechanical Engineering, University of Bath \
            
        ]
    }
)

#show raw.where(block: true): block.with(fill: luma(240), inset: 1em, radius: 0.5em, width: 100%)
#show figure: set block(breakable: true)
#show link: underline

// MARK: EX 1
= Introduction

*Finite Element Method (FEM)* is a powerful numerical technique for solving equations over a discrete domain. 
The simulated system is split into small regions called *elements*, connected by *nodes* which represent discrete points in the domain, together making up a *mesh*.
Elements are evaluated using *basis functions* which approximate the solution within each element based on node values @finite-element-method.
This approach allows for practical solutions to problems that may be difficult or impossible to solve analytically. 
In addition to this, the size and shape of elements can be adjusted to improve accuracy or reduce computational cost, making FEM a powerful and flexible tool for modelling (@fem-example).

#figure(
    image("resources/fem-example.jpg", width: 90%),
    caption: [Finite Element Modelling of a Wrench under a Test Load Scenario @comsol], 
)  <fem-example>

This coursework focuses on the implementation and verification of a FEM solver for the transient diffusion-reaction equation, given by @numerical-advection-diffusion-reaction:

#math.equation(
  block: true,
  numbering: "(1)",
  $ (delta c) / (delta t) = D (delta ^2 c)/(delta x^2) + lambda c + f $
) <transient-diffusion-reaction>

Where:
- $c$ is the concentration level
- $D$ is the diffusion coefficient
- $lambda$ is the reaction rate
- $f$ is a source term

The transient diffusion-reaction equation models processes where substances diffuse through a medium while undergoing reactions or being influenced by boundary interactions. 
Examples of situations modelled by this equation include the transfer of heat through a material or (as explored in Part 3 of this report) the diffusion of a drug through biological tissue.

This coursework describes the development and validation of a FEM solver for the transient diffusion-reaction equation. 
To keep the scope manageable, the solver was implemented in 1D, using MATLAB as the scripting language @matlab.

= Part 1: Software Verification

== Background

A static FEM solver was implemented in a previous coursework for the steady-state diffusion-reaction equation. 
This solver was subsequently adapted to solve the transient form of the equation (@transient-diffusion-reaction).

For the initial case, the values of $D = 1$ and $lambda = 0$ were used, representing a pure diffusion scenario with linear behaviour.
The *Crank-Nicolson* finite difference method was used for time integration. It has unconditional stability but no damping of oscillations, providing a good compromise between accuracy and stability at this stage @crank-nicolson.

The problem space was further defined with the following conditions:

#figure(
    caption: "Initial Case Conditions",
    block(width: 100%, inset: (top: 0%, bottom: 0%),
        align(center, //Align starts here
            table(
                columns: (auto, auto),
                inset: 5pt,
                align: horizon + left,
                [Problem Space], [$0 <= x <= 1$],
                [Left Boundary Condition], [Dirichlet: $c(0, t) = 0$],
                [Right Boundary Condition], [Dirichlet: $c(1, t) = 1$],
                [Initial Condition], [$c(x, 0) = 0$],
            )
        )
    )
) <intial-case-conditions>

These conditions have a known analytical solution, given by @analytical-solution-equation:

#math.equation(
  block: true,
  numbering: "(1)",
  $ c(x, t) = x + 2/pi sum^infinity_(n = 1) ((-1)^n)/n e^(-n^2 pi^2 t) sin(n pi x) $
) <analytical-solution-equation>

The analytical solution allows for direct comparison of results between the FEM solver and expected values, providing a quantitative measure of accuracy.

== Software Architecture

The solver was implemented with a modular, object-oriented software architecture to improve readability and control flow. 
Classes were created to encapsulate well-defined functions of the solver, such as mesh generation or plotting (@part1-software-architecture).

#figure(
    image("resources/part1/SoftwareArchitecture.drawio.png", width: 100%),
    caption: [High-Level Software Architecture of the FEM Solver], 
)  <part1-software-architecture>

== Results

Having implemented the FEM solver as described above, a simulation was run using a mesh size of 50 elements and a time step of 0.01s, over the time period $0 < t <= 1s$.

After this, the results were plotted on a series of charts for a visual comparison of the two solutions.
The first of these were heatmaps which are an effective method for visualising the 1D diffusion over time (@part1-analytical-heatmap, @part1-numeric-heatmap).

#figure(
    image("resources/part1/NumericHeatmap.png", width: 110%),
    caption: [FEM Solution of Diffusion Equation over using the Crank-Nicolson method over $0 <= x <= 1$ and \ $0 <= t <= 1s$], 
)  <part1-numeric-heatmap>

#figure(
    image("resources/part1/AnalyticalHeatmap.png", width: 110%),
    caption: [Analytical Solution of Diffusion Equation over $0 <= x <= 1$ and $0 <= t <= 1s$], 
)  <part1-analytical-heatmap>

The data was also represented in a 2D plot, showing the concentration through the mesh at sample times of $t = 0.05s, 0.1s, 0.3s, 1.0s$, shown in @part1-numeric-samples and @part1-analytical-samples.

Additionally, a chart was created for both solutions at a single point in the mesh ($x = 0.8$), shown in @analytical-x08. Unlike previous plots, this shows both methods on the same axes for direct comparison, demonstrating the agreement between the two solutions.

#figure(
    image("resources/part1/NumericSamples.png", width: 110%),
    caption: [FEM Solution of Diffusion Equation over using the Crank-Nicolson method over $0 <= x <= 1$ and at $t = 0.05s, 0.1s, 0.3s, 1.0s$],  
)  <part1-numeric-samples>

#figure(
    image("resources/part1/AnalyticalSamples.png", width: 110%),
    caption: [Analytical Solution of Diffusion Equation over $0 <= x <= 1$ and at $t = 0.05s, 0.1s, 0.3s, 1.0s$], 
)  <part1-analytical-samples>

#figure(
    image("resources/part1/BothX08.png", width: 110%),
    caption: [Comparison of Analytical and FEM Solutions at $x = 0.8$ over $0 <= t <= 1s$],  
)  <analytical-x08>

\

== Spacial and Temporal Convergence

To quantitatively assess the accuracy of the FEM solver, the *Root Mean Square (RMS)* error between numerical and analytical solutions was evaluated over a range of element and time step sizes.
As shown in @element-size-convergence and @element-temporal-convergence, the RMS error decreases with both smaller element sizes and smaller time steps, demonstrating convergence of the numerical solution towards the analytical solution with increasing resolution.

#figure(
    image("resources/part1/ElementSizeConvergence.png", width: 110%),
    caption: [Comparison of RMS errors at $t = 1s$ for Varying Element Sizes],  
)  <element-size-convergence>

#figure(
    image("resources/part1/TimeStepConvergence.png", width: 110%),
    caption: [Comparison of RMS errors at $t = 1s$ for Varying Time Steps],  
)  <element-temporal-convergence>

== Testing and Validation

A set of unit tests were created alongside the FEM solver, to verify the functionality of individual components such as mesh generation, element assembly, and time integration.
As part of the development process, the project was continuously tested to ensure it passed all scenarios.

In particular, a unit test was created to validate the solver against a manufactured solution of the transient diffusion-reaction equation.
This involved selecting specific values for $D$, $lambda$, and $f$ such that the solution could be expressed in a simple analytical form.


= Part 2: Software features

== Error Evaluation

In Part 1 of the coursework, the RMS error term was used to evaluate the accuracy of the FEM solver. 
While RMS is a useful metric, it can be sensitive to outliers and therefore may not always provide a complete picture of the solution accuracy.
L2 norm doesn't suffer as much from this, and is more widely used in literature as a result @numerical-advection-diffusion-reaction.
To address this, a dedicated L2 error evaluation class was added to the solver, allowing for more robust error analysis.

== Integration Methods

Using the L2 norm error evalutation class, the performance of three different time integration methods was compared: 
Forward (Explicit) Euler, Backward (Implicit) Euler, and Crank-Nicolson.

#figure(
    image("resources/part2/L2ErrorTimeIntegration.png", width: 110%),
    caption: [Comparison of RMS errors at $t = 1s$ for Varying Time Steps],  
)  <part2-time-integration-comparison>

After this, a study was run to compare the speed and stability of each method.


#figure(
    image("resources/part2/BasisComparison.png", width: 110%),
    caption: [Comparison of RMS errors at $t = 1s$ for Varying Time Steps],  
)  <part2-basis-function-comparison>


#figure(
    image("resources/part2/IntegrationComparison.png", width: 110%),
    caption: [Comparison of RMS errors at $t = 1s$ for Varying Time Steps],  
)  <part2-integration-comparison>


== Quadratic Basis Functions

== Gaussian Quadrature




- The *Mesh* and *MeshElement* classes were modified to support higher-order elements.

== Comparing Solver Methods

- Forward Euler - conditionally stable, requires 0.0002s dt for stability
   - needs to be dt <= (dx^2) / (2D)

- other two methods - unconditionally stable



\

Next, the FEM solver was extended to account for the following advanced features:
- Different solver methods (Explicit Euler, Implicit Euler, Crank-Nicolson)
- Gaussian quadrature for numerical integration
- Quadratic basis functions
- Using L2 norm to evaluate solution accuracy

\


NEED TO REFINE MESH!


= Part 3: Modelling & Simulation Results

= Conclusion

// MARK: REFERENCES
= References

#bibliography(
    "resources/ME40064-References.yml",
    title: none,
    style: "ieee"
)


= Use of Generative AI

This coursework was completed in Visual Studio Code (with the #link("https://marketplace.visualstudio.com/items?itemName=MathWorks.language-matlab", "MATLAB Extension")), using Typst for report writing. The #link("https://github.com/features/copilot", "GitHub Copilot") AI tool was enabled, providing generative suggestions for report phrasing and code snippets.

#pagebreak()

#set page(columns: 1)

= Appendix

#let sections = (
  ("Main", (
    "../src/main.m", 
    ""
  )),

  ("Coursework", (
    "../src/coursework/Coursework.m", 
    ""
  )),

  ("Mesh", (
    "../src/mesh/Mesh.m",
    "../src/mesh/MeshElement.m",
    "../src/mesh/MultilayerMesh.m",
    "../src/mesh/LayerProperties.m",
  )),

  ("Analytical", (
    "../src/analytical/AnalyticalSolver.m",
    "../src/analytical/TransientAnalyticSoln.m",
  )),

  ("Plotter", (
    "../src/plotter/Plotter.m", 
    ""
  )),

  ("Solution", (
    "../src/solution/Solution.m", 
    "../src/solution/L2Error.m", 
    "../src/solution/DoseEvaluator.m", 
    ""
  )),

  ("Solver", (
    "../src/solver/NumericSolver.m", 
    "../src/solver/BoundaryCondition.m", 
    "../src/solver/BoundaryType.m", 
    "../src/solver/ElementMatrices.m", 
    "../src/solver/IntegrationMethod.m", 
    "../src/solver/IntegrationType.m", 
    ""
  )),


  ("Tests", (
    "../src/tests/NumericSolverTest.m", 
    ""
  )),

)

#context [
    #let n = counter("code-count")

    #for (section, files) in sections {

        n.step()

        context [
            
            = #section

            #for file in files {
                if file != "" {
                    [
                        == #file.split("/").last()
                        #let text = read(file)
                        #raw(text, lang: "matlab", block: true)
                    ]
                }
            }

            #if n.get().at(0) < sections.len() [
                #pagebreak()
            ]
        ]
    }

]
