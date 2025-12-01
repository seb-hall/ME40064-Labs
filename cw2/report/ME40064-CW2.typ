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

This coursework focuses on the implementation and verification of a FEM solver for the transient diffusion-reaction equation, given by

#math.equation(
  block: true,
  numbering: "(1)",
  $ (delta c) / (delta t) = D (delta ^2 c)/(delta x^2) + lambda c + f, $
) <transient-diffusion-reaction>

where $c$ is the concentration level, $D$ is the diffusion coefficient, $lambda$ is the reaction rate
and $f$ is a source term @numerical-advection-diffusion-reaction.

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
The comparison  test was run using a mesh with 10 elements and a time step size of 0.0001s. 

#figure(
    image("resources/part2/L2ErrorTimeIntegration.png", width: 100%),
    caption: [Comparison of L2 Errors for Different Time Integration Methods],  
)  <part2-time-integration-comparison>

This shows that the Forward Euler method had a higher initial accuracy, approaching the solution more quickly than the other two methods, but that it started to decrease in accuracy again afterwards.
This was likely caused by instability in the method, as it is only conditionally stable.

To illustrate this further, a stability analysis was performed for all three methods, using a larger mesh of 50 elements. The stability results are shown in @part2-time-integration-stability-analysis:

#figure(
    caption: "Integration Method Stability Comparison",
    block(width: 100%, inset: (top: 0%, bottom: 0%),
        align(center, //Align starts here
            table(
                columns: (auto, auto, auto, auto),
                inset: 5pt,
                align: horizon + left,
                table.header(
                    [*dt*], [*Forward Euler*], [*Backward Euler*], [*Crank-Nicolson*]
                ),
                [0.0001], [Stable], [Stable], [Stable],
                [0.001], [Unstable], [Stable], [Stable],
                [0.01], [Unstable], [Stable], [Stable],
                [0.1], [Unstable], [Stable], [Stable],
                [0.25], [Unstable], [Stable], [Stable]
            )
        )
    )
) <part2-time-integration-stability-analysis>

This shows that the Forward Euler method was only stable for very small time steps, while the other two methods demonstrated *unconditional stability*, remaining stable across all tested time steps.

For linear finite elements, the stability condition for the Forward Euler method is given by the following equation @euler-stability:

#math.equation(
  block: true,
  numbering: "(1)",
  $ d t <= (d x^2) / (2D) $
) <euler-stability-equation>

Therefore, for a mesh with 50 elements over the domain $0 <= x <= 1$ and $D = 1$, the value of $d t$ must be no more than 0.0002s for stability, which aligns with the results shown in @part2-time-integration-stability-analysis.

== Gaussian Quadrature

So far, the solver has only been used with a simple 2-point trapezoidal integration method for evaluating element matrices.
While this method is easy to implement, it treats all elements as linear, requiring meshes with high numbers of elements to achieve good accuracy for non-linear problems.

*Gaussian Quadrature* is an alternative integration method that can provide a more accurate result with the same number of integration points as trapezoidal integration, resulting in a more efficient solution @gaussian-quadrature.

== Quadratic Basis Functions

For 2-point basis functions like those used in the coursework so far, Gaussian Quadrature with 2 points will produce an identical result to trapezoidal integration. 
The mesh was therefore updated to support higher-order basis functions, such as quadratic (3-point) elements, where each element has a node at each end and one in the middle.

The L2 error of a quadratic mesh with both trapezoidal and Gaussian integration methods is shown below in @part2-gaussian-trapezoidal-comparison:

#figure(
    image("resources/part2/L2ErrorGaussianTrapezoidal.png", width: 110%),
    caption: [Comparison of L2 Errors for Gaussian Quadrature and Trapezoidal Integration],  
)  <part2-gaussian-trapezoidal-comparison>

This shows a clear improvement in accuracy when using Gaussian Quadrature over trapezoidal integration with quadratic basis functions, approaching the analytical solution in a shorter time.

The reason for this improved performance is that Gaussian Quadrature evaluates the integrand at specially chosen points, capturing a more accurate representation of the function being integrated. @part2-gaussian-trapezoidal-diagram shows a visual comparison of a trapezoidal integration, alongside a Gaussian Quadrature with 2 intermediary points.

#figure(
    image("resources/part2/Comparison_Gaussquad_trapezoidal.svg.png", width: 100%),
    caption: [Visualisation of Gaussian Quadrature vs Trapezoidal Integration @gaussian-diagram],  
)  <part2-gaussian-trapezoidal-diagram>

In the code implementation, the mesh was implemented with an arbitrary order parameter, dynamically calculating the number of nodes based off the baseline element count and order. In comparison, the code for Gaussian Quadrature was implemented with pre-defined functions for 1, 2 and 3 integration points, and shape functions for linear and quadratic elements.

== Summary of Features

The addition of L2 error evaluation was an effective way to quantitatively assess the accuracy of the FEM solver, with varying configurations. It was found that the Crank-Nicolson method remained a suitable choice for time integration, balancing accuracy and stability, while the addition of Gaussian Quadrature and higher-order basis functions showed a significant improvement to solution accuracy.

Together these features enhance the capability and robustness of the FEM solver, allowing it to tackle a wider range of problems with improved accuracy and efficiency.

The robust, object-oriented software architecture introduced in Part 1 was also extended, with new classes and tests for the additional functionality.


= Part 3: Modelling & Simulation Results

== Overview

The transient FEM solver developed in Parts 1 and 2 was then applied to a practical problem: modelling the diffusion of a drug through a multilayer skin structure, as shown in the diagram below:

#figure(
    image("resources/part3/description-skin.png", width: 110%),
    caption: [1D Multilayer Finite Element Mesh of Skin Tissue Layers  @part3-description],  
)  <part3-skin-description>

\

The concentration of the drug is modelled by the transient diffusion-reaction equation

#math.equation(
  block: true,
  numbering: "(1)",
  $ (delta c) / (delta t) = D (delta ^2 c)/(delta x^2) - beta c - gamma c,$
) <part3-diffusion-reaction-equation>

where $c$ is the drug concentration, $D$ is the diffusion coefficient, $beta$ is the extra-vascular diffusivity, and $gamma$ is the drug degradation rate. For the purposes of modelling, $beta$ and $gamma$ are combined into a single reaction rate term i.e $lambda = beta + gamma$, as they both act as sink terms that reduce the drug concentration.

== Solver Modification

The main difference between the skin application and previous problems is the use of a multilayer mesh. A new *MultilayerMesh* class was created, inheriting from the original *Mesh* class, and overriding a method to generate a mesh made up of discrete layers (*MeshLayer* class), each with different properties. 

In addition to variable diffusion and reaction rates, a 'density ratio' property was added to each layer, allowing for a non-uniform distribution of elements across the mesh. Thinner layers can therefore be assigned a higher density ratio, resulting in a local mesh with higher resolution, and improved solution accuracy.

This was implemented in three passes. First, the total density of all layers was calculated. Then, the number of elements in each layer was found by multiplying the total element count by the ratio of the layer density to total density. As an example, if there were two layers with density ratios of 1 and 3 respectively, and a total of 40 elements, the layers would be assigned 10 and 30 elements respectively. After this, the node co-ordinates were generated for each layer in sequence, with a uniform distribution according to the number of elements assigned to that layer, and it's range of $x$ values.

== Simulation Results

The modified solver was then configured to use solve the coursework-specified problem, with the following conditions:

#figure(
    caption: "Drug Concentration Problem Conditions",
    block(width: 100%, inset: (top: 0%, bottom: 0%),
        align(center,
            table(
                columns: (auto, auto),
                inset: 5pt,
                align: horizon + left,
                [Problem Space], [$0 <= x <= 0.01$],
                [Left Boundary Condition], [Dirichlet: $c(0, t) = 30$],
                [Right Boundary Condition], [Dirichlet: $c(0.01, t) = 0$],
                [Initial Condition], [$c(x, 0) = 0$],
            )
        )
    )
) <part3-conditions>

This used a multilayer mesh with three layers representing the epidermis, dermis and sub-cutaneous tissue, with the following parameters:

#figure(
    caption: "Mesh Layer Parameters",
    block(width: 100%, inset: (top: 0%, bottom: 0%),
        align(center, 
            table(
                columns: (auto, auto, auto, auto),
                inset: 5pt,
                align: horizon + left,
                table.header(
                    [*Parameter*], [*Epidermis*], [*Dermis*], [*Sub-Cutaneous*]
                ),
                [$x$ Range], [$0 <= x < 0.00166667$], [$0.00166667 \ <= x < 0.005$], [$0.005 \ <= x <= 0.01$],
                [D], [4e-6], [5e-6], [2e-6],
                [$beta$], [0.0], [0.01], [0.01],
                [$gamma$], [0.02], [0.02], [0.02],
                [Density \ Ratio], [2.0], [1.0], [1.0],
            )
        )
    )
) <part3-mesh-parameters>

The simulation was then run using an initial mesh size of 50 elements and a time step of 0.01s, over the specified time period of $0 < t <= 30s$.

The results were plotted as a heatmap (@part3-initial-numeric-heatmap), showing the diffusion of the drug through the multilayer skin structure over time. Additionally marked on this plot are the approximate boundaries between each layer.

#figure(
    image("resources/part3/InitialNumericHeatmapMarkedup.png", width: 110%),
    caption: [FEM Solution of Drug Diffusion through \ Multilayer Skin Structure over $0 <= x <= 0.01$ and \ $0 <= t <= 30s$],  
)  <part3-initial-numeric-heatmap>

\

A stable profile is visible after around 10 seconds, with the epidermis layer almost immediately saturated to a high level, and the dermis soon after with a slightly lower concentration. The sub-cutaneous layer shows a much more gradual concentration gradient, mainly due to the Dirichlet boundary at $x = 0.01$ of $c(0.01, t) = 0$ forcing a perfect sink along the far edge of the mesh.

== Dose Evaluation

After establishing a baseline simulation profile, the next step was to calculate the effectiveness of the drug diffusion. This was evaluated using the *kappa* metric, defined as the integral of the concentration above a threshold level over a specified time period, given by 

#math.equation(
  block: true,
  numbering: "(1)",
  $ K = integral^(t = 30)_(t_op("eff")) c space d t $
) <part3-kappa-equation>

where $t_op("eff")$ is the time at which the concentration first exceeds the threshold level (here 4.0), $t=30$ is the end of the simulation period, and $c$ is the drug concentration at a target point in the mesh.

A new class, *DoseEvaluator*, was created to provide this functionality, with a static method *EvaluateSolution* that returns the kappa value for a given solution dataset, target location and threshold concentration. The logic is quite simple; first the node closest to the target location is identified, and then the concentration values for that node are checked against the target threshold. If the concentration exceeds the threshold, the solution is integrated between the time this occurs and the end of the simulation, using a trapezoidal method.

== Minimum Dose Search

The *DoseEvaluator* class was then modified to add a new function, *FindMinimumDose*. The purpose of this was to indentify the minimum required dose of the drug that would achieve a sufficient value of kappa. 

This process was achieved with a *binary search* algorithm, which iteratevely narrows down the range of possible dose values, converging on the minimum effective dose. It relies on having a correct initial upper and lower bound for the dose, but provides an effective search method with relatively few iterations.

The minimum effective dose, $c_op("DOSE")$, was defined as the value at which the concentration at $x = 0.005$ exceeds a threshold of $K > 1000$. A search was run with an initial search range of $0 <= c <= 100$, using a tolerance of 0.1 for convergence. The results of this search are shown below in @part3-dose-binary-search:

#figure(
    image("resources/part3/MinDoseBinarySearch.png", width: 100%), 
    caption: [Binary Search for Minimum Effective Dose],  
)  <part3-dose-binary-search>

This shows that the relationship between dose and kappa is linear, and that the minimum effective dose was found to be approximately *59.18*.


== Dose Sensitivity Analysis

A sensitivity analysis was then performed on the model, investigating the impact of varying the diffusion coefficient $D$, extra-vascular diffusivity $beta$, and drug degradation rate $gamma$ on the concentration at $x = 0.005$ over time, and the resulting dose effectiveness K (Kappa).

=== Diffusion Coefficient

The diffusion coefficient $D$ was investigated by scaling the original values for each layer by factors of 0.5, 0.75, 1.0, 1.5 and 2.0. The results of this analysis are shown below in @part3-diffusion-sensitivity and @part3-diffusion-kappa, illustrating the effect of varying $D$ on concentration over time, and of the resultant dose effectiveness. These plots show a clear trend of increasing diffusion coefficient leading to higher concentrations at the target point, and therefore higher dose effectiveness K. 

This is expected, as a higher diffusion coefficient allows the drug to spread more rapidly through the tissue layers, reaching the target point in a shorter time, with less degradation.

#figure(
    image("resources/part3/DiffusionSensitivityAnalysis.png", width: 100%), 
    caption: [Concentration at $x=D$ for Varying Values of $D$],  
)  <part3-diffusion-sensitivity>

#figure(
    image("resources/part3/DiffusionKappa.png", width: 100%), 
    caption: [Dose Effectiveness for Varying Values of $D$],  
)  <part3-diffusion-kappa>


=== Extra-Vascular Diffusivity

The next parameter to be investigated was extra-vascular diffusivity ($beta$), varied in the same way as the diffusion coefficient. The results of this analysis are shown below in @part3-beta-sensitivity and @part3-beta-kappa.

#figure(
    image("resources/part3/BetaSensitivityAnalysis.png", width: 100%), 
    caption: [Concentration at $x=D$ for Varying Values of $beta$],  
)  <part3-beta-sensitivity>

#figure(
    image("resources/part3/BetaKappa.png", width: 100%), 
    caption: [Dose Effectiveness for Varying Values of $beta$],  
)  <part3-beta-kappa>

These plots show an inverse, linear relationship between $beta$ and dose effectiveness. As $beta$ increases, the faster the drug diffuses out of the vascular system into surrounding tissue, reducing the concentration at the target point and therefore lowering K.

=== Drug Degradation Rate

Finally, the drug degredation rate ($gamma$) was investigated, again varied in the same way as previous parameters. The results of this analysis are shown below @part3-gamma-sensitivity and @part3-gamma-kappa.

#figure(
    image("resources/part3/GammaSensitivityAnalysis.png", width: 110%), 
    caption: [Concentration at $x=D$ for Varying Values of $gamma$],  
)  <part3-gamma-sensitivity>

#figure(
    image("resources/part3/GammaKappa.png", width: 110%), 
    caption: [Dose Effectiveness for Varying Values of $gamma$],  
)  <part3-gamma-kappa>

As with $beta$, these plots show an inverse, linear relationship between $gamma$ and dose effectiveness. A higher degradation rate results in the drug breaking down more quickly, reducing the concentration at the target point and lowering K. 

From a mathematical perspective, this is also expected as both $beta$ and $gamma$ act as sink terms in the diffusion-reaction equation, reducing the overall concentration. However, as the original values of $gamma$ were larger than those of $beta$, the impact of varying $gamma$ was more pronounced.

\

== Further Work

While the FEM solver developed in this coursework has proved effective for modelling the 1D drug diffusion problem, there are several areas where further work could improve it's performance.

One example is the implementation of continuous diffusion-reaction parameters across the mesh, instead of the discrete steps that are currently used. This would more accuractly represent the gradual changes in tissue properties that occur in real biological systems.

Alongside this, the boundary conditions of a perfect source at one side and a perfect sink at the other are idealised scenarios. More realistic boundary conditions, such as Neumann or Robin conditions that vary over time, could be implemented to better simulate real-world situations @neumann.

The values chosen for the model parameters were based on the values provided in the coursework brief. However, these were stated as approximate or in some places unrealistic. More researched and physically accurate values would therefore improve the validity of the simulation results.

Increasing the mesh or time fidelity would also be likely improve accuracy, although both result in higher computational cost. More advanced meshing techniques could be explored, taking the multi-density approach further by implementing adaptive techniques that refine the mesh further in areas of high gradient.

Finally, the solver could be extended to 2D or 3D problems @d1-d2, allowing for more complex geometries and diffusion scenarios to be modelled. This would require significant changes to the mesh generation and element assembly processes, but would greatly expand the range of applications for the solver.

= Conclusion

The FEM solver developed in this coursework was successfully implemented and validated against an analytical solution for the transient diffusion-reaction equation in Part 1. 
It was then improved in Part 2 with the addition of more advanced features such as L2 error evaluation, higher-order basis functions, and Gaussian Quadrature integration.

In Part 3, the solver was then applied to a practical problem, modelling the diffusion of a drug through a multilayer skin structure.
The initial results and subsequent analysis demonstrated a physically plausible diffusion profile, with a binary search method effectively identifying the minimum effective drug dose level.

Changes to key parameters in a sensitivity analysis produced results that aligned with both mathematical and physical interpretations of the sample problem, further validating the solver's performance.

Finally, several areas for further work were discussed and identified, providing a roadmap for future improvements to the solver. These further demonstrate the flexibility and suitability of the FEM approach for solving complex diffusion-reaction problems.

// MARK: REFERENCES
= References

#bibliography(
    "resources/ME40064-References.yml",
    title: none,
    style: "ieee"
)

#pagebreak()

#set page(columns: 1)

= Appendix - MATLAB Source Code

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
    "../src/mesh/MeshLayer.m",
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
