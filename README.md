# Hydrodynamic Billiards

The hydro_billiards repository contains code to perform simulations of electron hydrodynamics in 2D electron systems. The user can specify device geometry, effective temperature, and bulk scattering rates. The code performs a time-domain simulation of a fixed number of particles which potentially scatter off of eachother in energy/momentum conserving collisons. The success of a scattering event depends on phase-space considerations of the outgoing trajectories. Particularly it samples whether the outgoing trajectories would be into an occupied state.

In a real physical system, the outgoing states should be sampled from a fermi-dirac distrubtion, but for expedience here we set a minimum energy below which no simulated particle is allowed to occupy. Generally a combination of this energy minimum and the number of particles in the simulation map to the effective temperature, since temperature both sets the number of thermally activated quasiparticles as well as the available phase space to scatter into.

The example implementation of a hydrodome simulation was performed by navigating to the "example implmentation" folder in a terminal and executing:

`python ../runDomeScript.py dome800N_0p9E_10um.txt`

The parameters set by "dome800N_0p9E_10um.txt" are as follows:

- **Constriction Width** -- narrowest part of the injector (nominally in microns)

- **Scattering Probability** -- chance of direction being randomized in a given timestep. This maps to electron-phonon scattering

- **Minimum Energy** -- The minimum allowed energy a particle is allowed. By default particles start out averaging E = 1, so by setting a value lower than that, particle scattering is allowed

- **Source Drain Ratio** -- This is the fractional probability that an electron will be emitted by a source contact compared to a drain contact(probabilities are properly weighted based on the edge length). This is conceptually equivalent to applying a  voltage bias.

- **Number of Particles** -- fixed number of particles in a simulation. For convenience, a particle gets re-emitted as soon as it gets absorbed. This is equivalent to setting the quasi particle density, and is a component of setting the temperature.

- **time steps** -- number of iterations of simulation runs for. The default setting is to have particles propagate 10 nm per timestep.

- **save interval** -- the number of timesteps before the data is saved, overwriting the previous save. This is implemented to keep results in case of a program/system crash.

- **dome diameter** -- sets the diameter of the dome (nominally in microns)

- **injector height and width** -- sets parameters for the trapezoidal injector contact. (nominally in microns)

- **number of CPUs** -- sets the number of parallel simulations using the "multiprocessing" package. It runs this number of independant simulations. i.e. it gives more statistics, but not faster simulation time.

- **base output path** -- the path into which the output data will be stored



To generate an arbitrary geometry, follow the example in "example implementation/Jagged Rectangle". The geometry is defined by "jaggedRectangleCreator.ipynb" and the simulation is performed by navigating to the "example implmentation/Jagged Rectangle" folder in a terminal and executing:

`python ../../runBaseScript.py Jag_rectangle_4x20_700N_0p9E.txt`
