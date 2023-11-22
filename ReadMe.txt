
This directory accompanies that paper entitled "Accelerated Parallel Magnetic Resonance Imaging
with Compressed Sensing using Structured Sparsity" by Nicholas Dwork and Erin Englund.

To create results, there are two run files:
run_structuredSparseSENSE.m and
run_generateResults.m

run_structuredSparseSENSE simply generates a single image using the technique presented
in the manuscript with the wavelet transform.
run_generateResults creates an output directory and generates results for several reconstruction
algorithms for a set of regularization parameter values as well as a *.csv file with many
metrics of quality.  This code depends on CurveLab-2.1.3, which can be gotten from
www.curvelet.org.

Copyright, 2023

For inquiries, please contact Nicholas Dwork at nicholas.dwork@cuanschutz.edu


