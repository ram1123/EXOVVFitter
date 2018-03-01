<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [How to setup](#how-to-setup)
- [Some information](#some-information)

<!-- /TOC -->
# How to setup
To setup the code run following commands:

	cmsrel CMSSW_9_0_1
	cd CMSSW_9_0_1/src
	cmsenv
	git clone git@github.com:ram1123/EXOVVFitter.git
	cd EXOVVFitter
	scramv1 b clean; scramv1 b
	cd PDFs
	root -l compilePdfs.C


## Before running the code:

1. you must add the ntuples containing the same background (e.g. all the single top samples, as s-channel, t-channel, tW-channel.. into a single STop sample). <!-- In order to do that, you can use the hadda_mu.sh and hadda_el.sh scripts contained in the repository used to produce the small ntuples at the step2: -->
2. Provide the path of ntuples here: [link](https://github.com/ram1123/EXOVVFitter/blob/master/g1_exo_doFit_class_new.py#L146)
3. Change the root file read method if not running from store area of fnal from here: [link](https://github.com/ram1123/EXOVVFitter/blob/master/g1_exo_doFit_class_new.py#L2967)
3. If you modify one of the libraries inside the PDF/ folder, you must recompile it with the following commands: (fo instance, if you modify Util.cxx)

	cd PDFs/
	root -l compilePdfs.C


## To run the code:

	python g1_exo_doFit_class_new.py -b 
	

# Some information

1. `fit_AllSamples_Mj_and_Mlvj()` : This function fits signal and all background mj, mlvj (sideband) and mlvj (signal region) spectrum.
2. `get_data()` : Read data
3. `fit_WJetsNormalization_in_Mj_signal_region()` : Fits data's M\_j spectrum in sideband_low and sideband_hi to extract the WJets normalization in signal region.
4. `fit_mlvj_in_Mj_sideband()` : Fit data's M_lvj spectrum in sideband region and calculate the alpha from Wjets MC fitting result. So, that we can get the Wjets M_lvj shape in signal-region.
5. `get_WJets_mlvj_correction_sb_lo_to_signal_region()` :
6. `read_workspace()` : Draw the final M_lvj plots using the data_obs, signal and backgrounds shape and rate what saved in the root file.
7. `fit_mlvj_model_single_MC` :
