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
            *ME40064 System Modelling and Simulation - Coursework 3*
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
            *ME40064 System Modelling and Simulation - Coursework 3* \
        ]

        text(13pt)[
            #total-words Words, Candidate No. 11973, 9th January 2026\
            Department of Mechanical Engineering, University of Bath \
            
        ]
    }
)

#show raw.where(block: true): block.with(fill: luma(240), inset: 1em, radius: 0.5em, width: 100%)
#show figure: set block(breakable: true)
#show link: underline

// MARK: EX 1
= Introduction

Simulating dynamic systems can provide valuable insights into their behaviour under various conditions, but can also introduce challenges such as model accuracy and computation efficiency. 
For simulating vehicle dynamics, a *half car model*  (@halfcar) is often used to isolate the essential characteristics of a car's suspension system while simplifying the complexity of a full car model @half-car.
It models dynamics along a straight road, providing insight into vertical and pitch motions.

#figure(
    image("resources/halfcar.png", width: 100%),
    caption: [Representation of the Half Car Model @assignment-2-description],  
)  <halfcar>

The equations of motion for the half car model are as follows @assignment-2-description:

#math.equation(
  block: true,
  numbering: "(1)",
  $ m_B dot.double(s)_B =  F_text("front") + F_text("rear") - m_B g  $
) <halfcar-eqn-1>

#math.equation(
  block: true,
  numbering: "(1)",
  $ I dot.double(theta) = a F_text("front") - b F_text("rear")  $
) <halfcar-eqn-2>

Where $s_B$ is the vertical displacement, $theta$ is the pitch angle, $m_B$ is the body mass, $I$ is the moment of inertia, $g$ is gravitational acceleration, and $a$ and $b$ are distances from the center of mass to the front and rear axles respectively. 
$F_text("front")$ and $F_text("rear")$ are the forces exerted by the front and rear suspension models, described as a sprung mass in series with a damper and spring.

*Simulink* is a simulation environment and graphical modelling system, available as part of the MATLAB software suite @simulink.
It is particularly well-suited for modelling dynamic systems, providing a wide range of solvers and fundamental 'blocks' for constructing models (e.g. @mass-damper).

#figure(
    image("resources/massdamper.png", width: 100%),
    caption: [Simulink Spring-Mass-Damper System],  
)  <mass-damper>


= Part 1: Model Construction and Verification

== Creation of a Half Car Body Block

A Simulink block was created to solve the half car equations of motion. This had inputs of the front and rear suspension forces, and outputs of the body vertical displacement and pitch angle, as well as front and rear axle vertical displacements and velocities (@half-car-body).

#figure(
    image("resources/half-car-body-block.png", width: 90%),
    caption: [Simulink Half Car Body Block],   
)  <half-car-body>

The front and rear vertical displacements and velocities were calculated according to the following equations:

#math.equation(
  block: true,
  numbering: "(1)",
  $ s_text("front") = s_B + a theta, quad dot(s)_text("front") = dot(s)_B + a dot(theta) $
) <halfcar-eqn-3>

#math.equation(
  block: true,
  numbering: "(1)",
  $ s_text("rear") = s_B - b theta, quad dot(s)_text("rear") = dot(s)_B - b dot(theta) $
) <halfcar-eqn-4>

== Verification of the Half Car Body Block

A significant benefit of Simulink is the clear definition of inputs and outputs for each block, allowing for modular testing and verification of individual components.
For the half-car block, the only inputs were the front and rear suspension forces, alongside the set of parameters defined in the block mask.

The first set of tests focused on @halfcar-eqn-1, analysing the vertical displacement output, while the second set of tests focused on @halfcar-eqn-2, studying the pitch angle output.

=== Vertical Displacement Tests

First, the effect of the sum of the front and rear suspension forces on the vertical displacement output was investigated.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/car-body-tests/equal-weight.png", width: 100%)
            ]),
            caption: [$F_text("front")$ and $F_text("rear")$ Sum to $m_B g$],
        ) <fig-equal-weight>
    ],
    [

        #figure(
            block([
                #image("resources/car-body-tests/less-than-weight.png", width: 100%)
            ]),
            caption: [$F_text("front")$ and $F_text("rear")$ Sum to Less than $m_B g$],

        ) <fig-less-than-weight>
    ]
    )
])

In @fig-equal-weight, the front and rear suspension forces sum to equal the weight of the car body ($m_B g$), resulting in a balanced system with no change in vertical displacement, as expected from @halfcar-eqn-1. When the suspension force is reduced so that the resultant force is less than the weight of the car body, @fig-less-than-weight shows that the car body drops downwards, again as expected.

Next, the effect of the magnitude of the mass $m_B$ on the vertical displacement output was investigated, with a constant resultant force greater than the weight of the car body.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/car-body-tests/lower-mass.png", width: 100%)
            ]),
            caption: [1N Resultant Upward Force with 100kg Mass],
        ) <fig-lower-mass>
    ],
    [

        #figure(
            block([
                #image("resources/car-body-tests/higher-mass.png", width: 100%)
            ]),
            caption: [1N Resultant Upward Force with 1000kg Mass],

        ) <fig-higher-mass>
    ]
    )
])

@fig-lower-mass and @fig-higher-mass show that a lower mass results in a greater vertical displacement response to the same applied force, as expected from @halfcar-eqn-1. 

The distance the car body moves in response to the resultant 1N force matches the theoretical values calculated using Newton's second law ($F = m a$) over a 10 second interval. 

In addition to this, the suspension force being greater than the weight of the car body results in an upward acceleration, validating the model further.

=== Pitch Angle Tests

First, the effect of $a$ and $b$ on the pitch angle output was investigated.
The ratio of these parameters determines the position of the car's centre of gravity (CG), which should affect the pitch angle response to equal front and rear forces.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/car-body-tests/same-force.png", width: 100%)
            ]),
            caption: [Pitch Angle With Equal Front and Rear Forces, $a$ and $b$ Equal],
        ) <fig-unity-cg>
    ],
    [

        #figure(
            block([
                #image("resources/car-body-tests/forward-cg.png", width: 100%)
            ]),
            caption: [Pitch Angle With Equal Front and Rear Forces, CG Shifted Forward],

        ) <fig-forward-cg>
    ]
    )
])

@fig-unity-cg shows that the pitch angle remains at zero when equal front and rear forces are applied with $a$ and $b$ equal, as expected.
When the CG is shifted forwards, while keeping the front and rear suspension forces equal, @fig-forward-cg shows that the car pitches downwards as expected from @halfcar-eqn-2.

Next, the effect of unequal front and rear forces on the pitch angle output was investigated, with $a$ and $b$ kept equal.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/car-body-tests/higher-front.png", width: 100%)
            ]),
            caption: [Pitch Angle with Higher Front Force, $a$ and $b$ Equal],

        ) <fig-higher-front>
    ],
    [

        #figure(
            block([
                #image("resources/car-body-tests/higher-rear.png", width: 100%)
            ]),
            caption: [Pitch Angle with Higher Rear Force, $a$ and $b$ Equal],

        ) <fig-higher-rear>
    ] 
    )
])

@fig-higher-front shows that when the front suspension force is higher than the rear, the car pitches upwards as expected from @halfcar-eqn-2.
Conversely, @fig-higher-rear shows that when the rear suspension force is higher than the front, the car pitches downwards, again as expected.

Finally, the effect of the moment of inertia $I$ on the pitch angle output was investigated, with $a$ and $b$ kept equal and a higher front suspension force applied.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/car-body-tests/low-inertia.png", width: 100%)
            ]),
            caption: [Pitch Angle with Higher Front Force, Lower Moment of Inertia],

        ) <fig-low-inertia>
    ],
    [

        #figure(
            block([
                #image("resources/car-body-tests/high-inertia.png", width: 100%)
            ]),
            caption: [Pitch Angle with Higher Front Force, Higher Moment of Inertia],

        ) <fig-high-inertia>
    ] 
    )
])

@fig-low-inertia and @fig-high-inertia show that a lower moment of inertia results in a higher pitch angle response to the same applied forces, as expected from @halfcar-eqn-2.

\

#figure(
    caption: "Half Car Body Block Test Results",
    block(width: 100%, inset: (top: 0%, bottom: 0%),
        align(center, //Align starts here
            table(
                columns: (auto, auto),
                inset: 5pt,
                align: horizon + left,
                table.header(
                    [*Test*], [*Result*]
                ),
                [Body doesn't move with upward force of mg], [PASS],
                [Body moves down with upward force $<$ mg], [PASS],
                [Lower mass results in greater displacement], [PASS],
                [Balanced forces with CG at midpoint results in zero pitch], [PASS],
                [CG forward results in nose-down pitch], [PASS],
                [Higher front force results in nose-up pitch], [PASS],
                [Higher rear force results in nose-down pitch], [PASS],
                [Higher moment of inertia results in smaller pitch angle], [PASS],
                
            )
        )
    )
)

== Creation of Full Half Car Model Block

Having validated the half car body block, a full half car model was assembled by coupling it with front and rear suspension and tyre blocks, as shown in @halfcar-subsystem.

#figure(
    image("resources/half-car-subsystem.png", width: 110%),
    caption: [Simulink Half Car Model],  
)  <halfcar-subsystem>

The parameters of the suspension, tyre and body models were all set according to the coursework specification @assignment-2-description.

Simulink Goto/From blocks and bus signals were used to manage signal routing and maintain model clarity.

== Verification of the Full Half Car Model

As with the half car body block, a suite of tests was performed to verify the full half car model.

First, the settle response of the car and the effect of mass on settle displacement were investigated.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/half-car-tests/cw-settle.png", width: 100%)
            ]),
            caption: [Settle Response with Coursework Mass],

        ) <fig-cw-settle>
    ],
    [

        #figure(
            block([
                #image("resources/half-car-tests/half-mass-settle.png", width: 100%)
            ]),
            caption: [Settle Response with Half Mass],

        ) <fig-half-mass-settle>
    ] 
    )
])

@fig-cw-settle and @fig-half-mass-settle show that the resultant vertical displacement of the car after settling is correctly affected by the mass of the car body.

Next, the effect of suspension damping on the settle response was investigated.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/half-car-tests/settle-less-damping.png", width: 100%)
            ]),
            caption: [Settle Response with Decreased Damping],

        ) <fig-settle-less-damping>
    ],
    [

        #figure(
            block([
                #image("resources/half-car-tests/settle-more-damping.png", width: 100%)
            ]),
            caption: [Settle Response with Increased Damping],

        ) <fig-settle-more-damping>
    ] 
    )
])

@fig-settle-less-damping and @fig-settle-more-damping show that the damping coefficient of the suspension affects the rate of settling, as expected.

After this, the effect of the CG position on the pitch angle during settling was investigated.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/half-car-tests/cw-settle-pitch.png", width: 100%)
            ]),
            caption: [Settle Pitch with Coursework CG Position],

        ) <fig-cw-settle-pitch>
    ],
    [

        #figure(
            block([
                #image("resources/half-car-tests/equal-settle-pitch.png", width: 100%)
            ]),
            caption: [Settle Pitch with CG at Midpoint],

        ) <fig-unity-settle-pitch>
    ] 
    )
])

@fig-cw-settle-pitch and @fig-unity-settle-pitch show that under the coursework parameters, the car settles with a small pitch angle due to the CG being forward of centre, with identical front and rear suspension characteristics. When the CG is moved to the midpoint between the front and rear axles, the car settles with zero pitch angle, as expected.

Next, the response of the car to road steps was investigated, first with a step hitting both front and rear axles simultaneously, and then with a staggered step hitting the front axle before the rear.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/half-car-tests/same-bump-s.png", width: 100%)
            ]),
            caption: [Displacement Response to Step],

        ) <fig-same-bump-s>
    ],
    [

        #figure(
            block([
                #image("resources/half-car-tests/same-bump-theta.png", width: 100%)
            ]),
            caption: [Pitch Response to Step],

        ) <fig-same-bump-theta>
    ],
    [
        #figure(
            block([
                #image("resources/half-car-tests/staggered-bump-s.png", width: 100%)
            ]),
            caption: [Displacement Response to Staggered Step],

        ) <fig-staggered-bump-s>
    ],
    [

        #figure(
            block([
                #image("resources/half-car-tests/staggered-bump-theta.png", width: 100%)
            ]),
            caption: [Pitch Response to Staggered Step],

        ) <fig-staggered-bump-theta>
    ] 
    )
])

These results show that the car responds appropriately to road disturbances.

#figure(
    caption: "Half Car Model Test Results",
    block(width: 100%, inset: (top: 0%, bottom: 0%),
        align(center, //Align starts here
            table(
                columns: (auto, auto),
                inset: 5pt,
                align: horizon + left,
                table.header(
                    [*Test*], [*Result*]
                ),
                [Settle displacement affected by mass], [PASS],
                [Lower damping results in slower settling], [PASS],
                [Higher damping results in faster settling], [PASS],
                [CG forward results in nose-down pitch during settling], [PASS],
                [CG at midpoint results in zero pitch during settling], [PASS],
                [Step input results in appropriate vertical and pitch response], [PASS],
                [Staggered step input results in appropriate vertical and pitch response], [PASS],
                
            )
        )
    )
)


= Part 2: Investigation of Car \ Performance

A sinusoidal road profile was created to investigate the performance of the half car model under different conditions.
The road profile had an amplitude of 0.01m and a wavelength of 1m, and was tested with a range of speeds. 

A baseline dataset was captured, for the core performance metrics of body vertical acceleration ($dot.double(s)_B$), body pitch angular acceleration ($dot.double(theta)$), and front wheel vertical displacement relative to the road ($s_(w, text("front"))$), at speeds of 1m/s and 5m/s.

== Baseline Performance

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/perf/sine/1ms-sddot.png", width: 100%)
            ]),
            caption: [$dot.double(s)_B$ at 1m/s, Nominal Parameters],

        ) <fig-1ms-sddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine/5ms-sddot.png", width: 100%)
            ]),
            caption: [$dot.double(s)_B$ at 5m/s, Nominal Parameters],

        ) <fig-5ms-sddot>
    ],
    [
        #figure(
            block([
                #image("resources/perf/sine/1ms-thetaddot.png", width: 100%)
            ]),
            caption: [$dot.double(theta)_B$ at 1m/s, Nominal Parameters],

        ) <fig-1ms-thetaddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine/5ms-thetaddot.png", width: 100%)
            ]),
            caption: [$dot.double(theta)_B$ at 5m/s, Nominal Parameters],

        ) <fig-5ms-thetaddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine/1ms-sw.png", width: 100%)
            ]),
            caption: [$s_(w, text("front"))$ at 1m/s, Nominal Parameters],

        ) <fig-1ms-sw>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine/5ms-sw.png", width: 100%)
            ]),
            caption: [$s_(w, text("front"))$ at 5m/s, Nominal Parameters],

        ) <fig-5ms-sw>
    ] 
    )
])

== Performance with Increased \ Suspension Stiffness

After establishing the baseline performance, the suspension stiffness was doubled to investigate its effect on the performance metrics.

As specified by the coursework description @assignment-2-description, the spring contants of both front and rear suspensions used a lookup table, so all values in the table were multiplied by 2 to achieve the increased stiffness.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/perf/sine-2ks/1ms-sddot.png", width: 100%)
            ]),
            caption: [$dot.double(s)_B$ at 1m/s, $2k_s$],

        ) <fig-2ks-1ms-sddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2ks/5ms-sddot.png", width: 100%)
            ]),
            caption: [$dot.double(s)_B$ at 5m/s, $2k_s$],

        ) <fig-2ks-5ms-sddot>
    ],
    [
        #figure(
            block([
                #image("resources/perf/sine-2ks/1ms-thetaddot.png", width: 100%)
            ]),
            caption: [$dot.double(theta)_B$ at 1m/s, $2k_s$],

        ) <fig-2ks-1ms-thetaddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2ks/5ms-thetaddot.png", width: 100%)
            ]),
            caption: [$dot.double(theta)_B$ at 5m/s, $2k_s$],

        ) <fig-2ks-5ms-thetaddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2ks/1ms-sw.png", width: 100%)
            ]),
            caption: [$s_(w, text("front"))$ at 1m/s, $2k_s$],

        ) <fig-2ks-1ms-sw>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2ks/5ms-sw.png", width: 100%)
            ]),
            caption: [$s_(w, text("front"))$ at 5m/s, $2k_s$],

        ) <fig-2ks-5ms-sw>
    ] 
    )
])

The stiffer suspension resulted in a smoother ride at lower speeds, with reduced body accelerations (as shown in @fig-1ms-sddot vs @fig-2ks-1ms-sddot, and @fig-1ms-thetaddot vs @fig-2ks-1ms-thetaddot), alongside lower wheel displacement relative to the road (@fig-1ms-sw vs @fig-2ks-1ms-sw).

At higher speeds, however, the stiffer suspension led to increased body accelerations (@fig-5ms-sddot vs @fig-2ks-5ms-sddot, and @fig-5ms-thetaddot vs @fig-2ks-5ms-thetaddot), indicating a harsher ride. The wheel displacement relative to the road also increased (@fig-5ms-sw vs @fig-2ks-5ms-sw), suggesting reduced handling characteristics.

The improved performance at lower speeds can be attributed to the stiffer suspension's ability to better resist body movement in response to changes in road profile, while the degraded performance at higher speeds likely results from a shift in the suspension's natural frequency, leading to an increased amplitude compared to the baseline.

== Performance with Increased \ Suspension Damping

Next, the suspension damping was doubled to investigate its effect on the performance metrics.
Also specified by the coursework description @assignment-2-description, the damping coefficients were different depending on if the suspension was in compression or extension.

Both the compression and extension damping coefficients of the front and rear suspensions were therefore multiplied by 2 to achieve the increased damping for analysis.

#block(width: 100%,
inset: (left: 0%, right: 0%),
[
    #grid(
    columns: (0.5fr, 0.5fr),
    gutter: 10pt,
    align: horizon,
    [
        #figure(
            block([
                #image("resources/perf/sine-2cs/1ms-sddot.png", width: 100%)
            ]),
            caption: [$dot.double(s)_B$ at 1m/s, $2c_s$],

        ) <fig-2cs-1ms-sddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2cs/5ms-sddot.png", width: 100%)
            ]),
            caption: [$dot.double(s)_B$ at 5m/s, $2c_s$],

        ) <fig-2cs-5ms-sddot>
    ],
    [
        #figure(
            block([
                #image("resources/perf/sine-2cs/1ms-thetaddot.png", width: 100%)
            ]),
            caption: [$dot.double(theta)_B$ at 1m/s, $2c_s$],

        ) <fig-2cs-1ms-thetaddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2cs/5ms-thetaddot.png", width: 100%)
            ]),
            caption: [$dot.double(theta)_B$ at 5m/s, $2c_s$],

        ) <fig-2cs-5ms-thetaddot>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2cs/1ms-sw.png", width: 100%)
            ]),
            caption: [$s_(w, text("front"))$ at 1m/s, $2c_s$],

        ) <fig-2cs-1ms-sw>
    ],
    [

        #figure(
            block([
                #image("resources/perf/sine-2cs/5ms-sw.png", width: 100%)
            ]),
            caption: [$s_(w, text("front"))$ at 5m/s, $2c_s$],

        ) <fig-2cs-5ms-sw>
    ] 
    )
])

As with the increased stiffness case, the increased damping resulted in a smoother ride at lower speeds, with reduced body accelerations (as shown in @fig-1ms-sddot vs @fig-2cs-1ms-sddot, and @fig-1ms-thetaddot vs @fig-2cs-1ms-thetaddot), alongside lower wheel displacement relative to the road (@fig-1ms-sw vs @fig-2cs-1ms-sw).

Also like the increased stiffness case, at higher speeds the increased damping led to increased body accelerations (@fig-5ms-sddot vs @fig-2cs-5ms-sddot, and @fig-5ms-thetaddot vs @fig-2cs-5ms-thetaddot), indicating a harsher ride. The wheel displacement relative to the road also increased (@fig-5ms-sw vs @fig-2cs-5ms-sw), and my more than in the increased stiffness case, suggesting further reduced handling characteristics.

At lower speeds, the improved performance is likely due to the increased energy dissipation provided by the higher damping, reducing oscillations in response to road profile changes. At high speeds, however, the dampers may be over-damping the system, making it act more like a rigid body and transmitting more road disturbances to the sprung mass.

= Further Work

= Conclusion

Talk about adaptive dampers?


// MARK: REFERENCES
= References

#bibliography(
    "resources/ME40064-References.yml",
    title: none,
    style: "ieee"
)

#pagebreak()

#set page(columns: 1)

= Appendix - Simulink Source Code
