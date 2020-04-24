# CATS
CATS clustering algorithm source codes and howto guide

Jacob Ezerski / Cheung Group
University of Houston

Introduction: 
Combinatorial Averaged Transient Structure (CATS) is a clustering algorithm designed for use with proteins that experience a high degree of dynamic behavior. A detailed description of the algorithm and case-study results can be found in the original publication: https://doi.org/10.1021/acs.jpcb.8b08852 
The clustering process requires several steps to reduce the original trajectory data into clusters. As with other clustering algorithms, a clustering coordinate must be used (RMSD, Rg, Etc..). CATS differs from many standard algorithms by using coordinates that display Gaussian distributions. Any coordinate set that contains Gaussian-like distributions may be used for clustering in CATS, however the codes included here and the rest of this manual was formulated using protein dihedral Phi and Psi angles.

Procedure:

1. Reduction of MD trajectory to dihedral coordinates: The first step is a generation of a dihedral coordinate data file from the raw trajectory. This can be done several different ways, however the output format must contain the phi and psi dihedral coordinates in columns for each residue in the target protein. Each line of the coordinate file must correspond to a frame from the trajectory. Ex: If you wish to analyze a 20-residue protein using a trajectory of 1000 frames, the data file will contain 1000 lines with 40 columns (separated by spaces). We have included a TCL script for generating this file using VMD named "getdihedrals.tcl"
With a trajectory loaded into VMD, open the TCL console and enter >source getdihedrals.tcl. The output of this program will be the dihedral coordinate file "dihedralangles.out"

2. Histogram of dihedral angles: A histogram of the dihedral angle coordinates can be generated using the "makeHistogram" program. It can be compiled by running the following line 
g++ makeHistogram.cpp -o makeHistogram.exe
The input of this program is the dihedral coordinate file generated from step 1. makeHistogram will prompt the user for the number of coordinates being analyzed. For a protein with N residues, the user should enter 2N at the prompt and press enter. The program will output a comma-separted histogram of the dihedral coordinates using 100 bins of size 3.6 degrees (since dihedrals are periodic). The histogram will contain 2N columns (where N is the number of residues) and 100 lines corresponding to the histogram bins. This output is named "histogramOutput.csv"

3. Gaussian curve fitting: The each coordinate distribution from the histogram generated in step 2 must be fitted with a Gaussian curve. Short trajectories or proteins with no stable structures will produce poor distributions. A matlab code "catsinputgenerator.m" can be used for fitting the trajectory with Gaussian curves and producing the probabilistic file outputs needed for the next step. In matlab, import the histogram file from step 2 and run catsinputgenerator.m with the histogram as input. For each coordinate, the user will be propted to identify the number of Gaussian-like peaks, choose a starting point for curve fitting, and accept the approximate Gaussian fit. When choosing a starting point, use the mouse to click at a local minima on the graph. Be sure that all of the fitted curves are correct otherwise clustering may not work. Once all coordinate distributions have been analyzed, there will be two output files: CATS_probabilities.out and CATS_deviations.out
These files contain the probabilities, averages and standard deviations associated with each coordinates dihedral distribution. They are used as inputs in the next step.

4. Clustering: The program "frame_extract" performs the final analysis. input arguments are: frame_extract [traj_dihedrals] [deviation file] [probfile] [epsilon] [#of top clusters] [% (optional)]

The first argument uses the dihedral data file from step 1. Argument 2 and 3 use the CATS_deviations.out and CATS_probabilities.out files generated in step 3. CATS compares the dihedral coordinates in the trajectory to the Gaussian curves generated in step 3 to determine what trajectory frames/conformations are similar. Epsilon multiplies the standard deviation of each coordinate distribution -- epsilon=1 will result in more clusters with each cluster having a smaller population, however the clusters will contain closer matched structures, while epsilon=3 will yield fewer clusters with larger cluster sizes and more structure deviation within each cluster. 
Enter the number of top clusters, which will output a list of the top M clusters (based on population size). Alternatively, enter the % to output a list of clusters accounting for X% of the trajectory. 

NOTES: Cluster relaxation is an option to allow clusters to form even though a small number of the coordinate distributions do not match perfectly. 

Coordinates can also be ignored in the analysis by adding # at the start of the line in the CATS_probabilities.out file. This is useful of the N or C-terminal residues fluctuate and cause many additional clusters to be generated. 

Output of the program will describe the top clusters, what frames in the trajectory the center structures correspond to, and the population of the cluster. Output files detail all clusters. 
