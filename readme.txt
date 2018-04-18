This GC compartmental model accompanies the manuscript 
Beining et al (2017): T2N as a new tool for robust electrophysiological modeling demonstrated for mature and adult-born dentate granule cells. eLife

Note! Be sure to add the Modules subfolder to your Matlab path, as this folder contains the T2N and treestoolbox packages!

T2N automatically compiles the .mod files in the lib_mech folder necessary for the model to be working. 

Important functions/scripts:
GC_experiments.m	is the main script that initializes the model, runs the experiments and makes figures


Scripts for understanding:
GC_exp_tests.m		is similar to "GC_experiments.m" but can be better used for playing with the model
GC_initModel.m		initializes the model after the parameters were defined in the first section of GC_experiments.m or GC_exp_tests.m
GC_biophys.m 		comprises the definition of all the channels and channel densities
GC_spinedensity.m	comprises the definition of spine densitiy scaling factors

Scripts/functions that make or plot experiments start with "aGC_" and are all incorporated in "GC_experiments.m" or at least "GC_exp_tests.m"

Notice that you have to recreate the dll file once you change something in an mod file.