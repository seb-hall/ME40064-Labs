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
            #total-words Words, Candidate No. 11973, 4th November 2025\
            Department of Mechanical Engineering, University of Bath \
            
        ]
    }
)

#show raw.where(block: true): block.with(fill: luma(240), inset: 1em, radius: 0.5em, width: 100%)
#show figure: set block(breakable: true)
#show link: underline

// MARK: EX 1
= Introduction

Finite Element Method (FEM) is a powerful numerical technique for solving partial differential equations (PDEs) over a discrete domain. This coursework focuses on the implementation and verification of a FEM solver for the transient diffusion-reaction equation, given by:

$ (delta c) / (delta t) = D (delta ^2 c)/(delta x^2) + lambda c + f $

= Part 1: Software Verification

== Overview

For this section, a static FEM solver was adapted to solve the transient diffusion equation (i.e #sym.lambda = 0 and f = 0), to give an equation of the form:

$ (delta c) / (delta t) = D (delta ^2 c)/(delta x^2) $

Furthermore, the domain was specified with:
- An x range of 0 to 1
- Dirichlet boundary conditions of c(0, t) = 1 and c(1, t) = 0
- An initial condition of c(x, 0) = 0

== Implementation

== Results

The solver was implemented and the following plots generated:


#figure(
    image("resources/AnalyticalHeatmap.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results], 
)  <analytical-heatmap>


#figure(
    image("resources/SolverHeatmap.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results], 
)  <solver-heatmap>

#figure(
    image("resources/AnalyticalSamples.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results], 
)  <analytical-samples>

#figure(
    image("resources/NumericSamples.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results],  
)  <solver-samples>


#figure(
    image("resources/AnalyticalX08.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results],  
)  <analytical-x08>


#figure(
    image("resources/NumericX08.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results],  
)  <solver-x08>

== Validation

= Part 2: Software features

Next, the FEM solver was extended to account for the following advanced features:
- Different solver methods (Explicit Euler, Implicit Euler, Crank-Nicolson)
- Gaussian quadrature for numerical integration
- Quadratic basis functions
- Using L2 norm to evaluate solution accuracy

\

- The *Mesh* and *MeshElement* classes were modified to support higher-order elements.

== Comparing Solver Methods

- Forward Euler - conditionally stable, requires 0.0002s dt for stability
   - needs to be dt <= (dx^2) / (2D)

- other two methods - unconditionally stable


\

NEED TO REFINE MESH!


= Part 3: Modelling & Simulation Results

= Conclusion

= Use of Generative AI

This coursework was completed in Visual Studio Code (with the #link("https://marketplace.visualstudio.com/items?itemName=MathWorks.language-matlab", "MATLAB Extension")), using Typst for report writing. The #link("https://github.com/features/copilot", "GitHub Copilot") AI tool was enabled for this, and provided generative suggestions for code snippets and report phrasing.
