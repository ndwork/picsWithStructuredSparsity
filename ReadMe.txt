
This directory accompanies that paper entitled "Accelerated Parallel Magnetic Resonance Imaging
with Compressed Sensing using Structured Sparsity" by Nicholas Dwork and Erin Englund.

To create results, there are two run files:
run_structuredSparseSENSE.m and
run_generateResults.m

The data for these codes to function is shared in the SPIRiT v0.3 directory by Michael Lustig
at https://people.eecs.berkeley.edu/~mlustig/Software.html.  Copy the brain_8ch.mat
file from that download into the data subdirectory, and then hit run.

run_structuredSparseSENSE generates a single image comparing the technique presented in the
manuscript with other existing techniques.
run_generateResults creates an output directory and generates results for several reconstruction
algorithms for a set of regularization parameter values as well as a *.csv file with many
metrics of quality.  This code depends on CurveLab-2.1.3, which can be gotten from
www.curvelet.org.

Copyright, 2023

For inquiries, please contact Nicholas Dwork at nicholas.dwork@cuanschutz.edu


