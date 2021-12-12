# ActiveDumbbellSim_Py

Code to simulate a system of Dumbbells in a fluid driven by a fixed magnitude force directed along its axis. 
Motivated by the following papers:

* Matan Yah Ben Zion, Alvin Modin, Yaelin Caba, Paul M Chaikin, *Cooperation in a fluid swarm of fuel-free micro-swimmers*, arXiv:2012.15087
* Antonio Suma, Giuseppe Gonnella, Gianluca Laghezza, Antonio Lamura, Alessandro Mossa, and Leticia F. Cugliandolo, *Dynamics of a homogeneous active dumbbell system*, Phys. Rev. E 90, 052130

## About the Implementation
The Code is implemented in Cython for performance. 

Truncated and Shifted Lennard-Jones Potential implemented for paired interactions.

Verlet Neighbour Lists with automatic update scheme implemented for efficient calculations of paired forces.

To incorporate effects of Brownian Motion of fluid particles on the Dumbbells and Viscosity,
the system is described by Langevin Equations.
Time Propagation Algorithm based on the following paper:

* Eric Vanden-Eijnden, Giovanni Ciccotti, *Second-order integrators for Langevin equations with holonomic constraints*, Chemical Physics Letters, Volume 429, Issues 1–3, 2006

In the absence of temperature/viscosity, equations are identical to Velocity Verlet.


## About the Code and Repo
As of now, the code is contained in Jupyter Notebooks for initial development and testing.
Each notebook is self-contained, corresponding to different implementations.
Eventually the code will moved out of the notebooks into respective modules.

CPU parallelization will be implemented soon.

Integrator, Pair Interactions and Neighbour Lists tested for accuracy using simulations of LJ Fluids. A more extensive benchmark will be performed soon after optimization/parallelization.

Also, a folder containing example simulations (png and mp4 files) and presentation pdf also in repo for reference. Will be replaced with proper documentation later. 
