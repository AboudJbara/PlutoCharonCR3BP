## Build Process

The goal of this project was to model and visualize the zero-velocity curves, equilibrium points, and orbital behaviors, such as tadpole and horseshoe orbits, in the Pluto-Charon system, as well as in a simplified, classical low-mass-ratio model. 

Everything was nondimensionalized. The primary separation, total mass, and angular rate were all set to 1. The primaries were placed at (-μ, 0) and (1 - μ, 0). The project uses two separate values of μ, the first is the real Pluto-Charon mass ratio (~0.11), and a low-mass ratio to show traditional orbits (0.001)

The potential and its gradients were defined with the function ```omegaValues()``` and the Jacobi constant was derived from the formula C = 2Ω - (vx² + vy²). The collinear Lagrange points were found with Newton-Raphson iteration, with each point having its own distinct branch equation, since each one requires a different sign convention. L4 and L5 were defined analytically, with known formulas. 

The equations of motion were integrated using a fourth-order Runge-Kutta (RK4) method. Earlier versions would drift due to step-size errors, but this was later fixed. The propagation loop records positions, checks if there is a collision, and creates arrays that represent the trajectory of the massless particle. The time step (h) was kept to 5 x 10⁻⁴ throughout scripts (although it was decreased to 5 x 10⁻³ in the master file for quicker processing times). The total duration ranged from 30 to 300, and can be varied to represent different types of orbits.

Most issues came from indexing and contouring. The Z-Grid initially mixed up i and j ordering, and flipped or displaced contours. Another early problem was plotting contours at L4 or L5, but the program would not display anything. This was fixed by adding a small offset to the C (1 x 10⁻³), since the contour at L4 and L5 is a minuscule point, and cannot be seen regularly. 

Choosing the right energy was critical for seeing specific behaviors. Horseshoes only appeared when C was set just below the C of L3, while tadpoles required energy near the C of L4. With the Pluto-Charon ratio, the same setup produced chaotic motion instead of stability. This is an accurate physical difference, and not an error.

A μ < 0.039 was required for a stable Horseshoe, and after various attempts to reproduce one at μ = 0.01, the trajectory would either stabilize into a tadpole or become chaotic. This was fixed after changing μ to 0.001. Both values of μ are under the critical ratio (0.039); however, 0.01 was insufficient, since the window for horseshoes was difficult to locate, so a smaller mass ratio was used for easier access to horseshoe orbits.

Each .py file in /scripts corresponds to a figure, with comments describing which small parameter changes reproduce the others. The standalone master file reproduces the main horseshoe case by default but can switch to Pluto-Charon with a single value change.

Validation was done by checking that the Jacobi drift stayed below 10⁻⁵, and that halving the step size produced consistent results. We also ensured that trajectories never crossed their own zero-velocity boundaries. The results for a small μ were smooth and conserved, while for the Pluto-Charon μ were irregular and chaotic, confirming the code's correctness and sensitivity.

This simulator is limited to the planar, circular case of the Pluto-Charon system. It ignores the eccentricity of orbits, although the structure can be extended to include them later. The final code provides a compact, self-contained implementation that reproduces textbook CR3BP behavior and demonstrates the transition from stability to chaos with an increasing mass ratio.
