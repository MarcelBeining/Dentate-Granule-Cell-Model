Before you can use the model, you have to compile the .mod files in the lib_mech folder into nrnmech.dll using the "mknrndll" function of NEURON (can be found in the Windows Start menue normally). Then, be sure that the name of the created dll-file is the same as defined in "params.nrnmech" in the first section of "GC_experiments.m" or "GC_exp_tests.m". 
Note! Many scripts from the model also need functions that are located in the "general code" folder in the parent directory (including T2N of course). Be sure that you added this folder to your Matlab path!

The important functions/scripts:
Model folder:
GC_experiments.m	is the main script that initializes the model, runs the experiments and makes figures

T2N folder:
t2n_writetrees.m	function needed to create the hoc file for a morphology. Once the hoc file is created and the mtr file is resaved
			(in order to save an ID leading to the hoc file, you will be asked for saving), the hoc file does not need to be created 
			anymore, unless you change something in the morphology.


Scripts for understanding:
Model Folder:
GC_exp_tests.m		is similar to "GC_experiments.m" but can be better used for playing with the model
GC_initModel.m		initializes the model after the parameters were defined in the first section of mGC_experiments.m or mGC_SH_test2.m
GC_biophys.m 		comprises the definition of all the channels and channel densities
GC_spinedensity.m	comprises the definition of spine densitiy scaling factors

Scripts/functions that make or plot experiments start with "aGC_" and should all be incorporated in "C_experiments.m" or at least "GC_exp_tests.m"


Tipps:
1) You may use the Matlab command "rename_nrnmech(X)" to automatically rename the nrnmech.dll to "nrnmech_winX.dll". This feature is useful if you work with the same folder on different machines that need different compiled nrnmech.dll files and you do not want to recreate the dll each time you change. Instead you can simply change the name in params.nrnmech, e.g. from "nrnmech_win1.dll" into "nrnmech_win2.dll" once you created and renamed both dll files.
2) Notice that you have to recreate the dll file(s) once you change something in an mod file or add a new one.



