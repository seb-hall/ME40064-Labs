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
            #total-words Words, \Candidate No. 11973, 4th November 2025\
            Department of Mechanical Engineering, University of Bath \
            
        ]
    }
)

#show raw.where(block: true): block.with(fill: luma(240), inset: 1em, radius: 0.5em, width: 100%)
#show figure: set block(breakable: true)
#show link: underline

// MARK: EX 1
= Introduction

= Part 1: Software Verification
_Extend FEA code to solve the transient form of the diffusion-reaction equation._

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
    image("resources/SolverSamples.png", width: 110%),
    caption: [Linear Reaction Operator Unit Test Results],  
)  <solver-samples>


= Part 2: Software features

= Part 3: Modelling & Simulation Results

= Conclusion

= Use of Generative AI

This coursework was completed in Visual Studio Code (with the #link("https://marketplace.visualstudio.com/items?itemName=MathWorks.language-matlab", "MATLAB Extension")), using Typst for report writing. The #link("https://github.com/features/copilot", "GitHub Copilot") AI tool was enabled for this, and provided generative suggestions for code snippets and report phrasing.
