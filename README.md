# btp-granular-rot-cyl
A repository of codes and files used for my B. Tech project (under the guidance of Prof. Devang Khakhar) in my final semester at IIT Bombay. The problem statement dealt with modelling the dynamics of granular particles in a partly-filled rotating cylinder.

The basic workflow was as follows:
 - LAMMPS simulations
 - Visualisation in Ovito
 - Analysis and post-processing using Python

Primarily to save myself the effort of running expensive computations again, the ```codes``` folder already contains the simulation output (the ```v.dat``` file). This enables the user to directly use the ```velocity.ipynb``` file for all the analysis. The ```rc.dump``` file to visualise the system in Ovito is unfortunately too large to upload, and will require re-running the simulations using the bash scripts.
