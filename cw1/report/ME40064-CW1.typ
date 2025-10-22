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

// MARK: INTRO
= Introduction

// MARK: EX 1
= Q1 - Local Element Matrix Functions

== Part A: MATLAB 2x2 Element Matrix


=== Code
```matlab
function matrix = LaplaceElemMatrix(D, eID, msh)
%% LaplaceElemMatrix - calculates a 2x2 element matrix for the diffusion 
%% operator in a 1d finite element mesh.
    
    % create base matrix
    matrix = [1, -1; -1, 1];    

    % calulate element size
    elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1);

    % scale base matrix by element size
    matrix = matrix * (D / elemSize);
 
end
```

=== Test Result
#figure(
    image("resources/q1a_test_result.png", width: 100%),
    caption: [Unit test results],
)  <and-network-diagram>



// MARK: CONCLUSION
== Possible Improvements


// MARK: REFERENCES
= References

#bibliography(
    "resources/ME40064-References.yml",
    title: none,
    style: "ieee"
)