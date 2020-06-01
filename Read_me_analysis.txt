----------------------------------------------------------------
Date created: 1st june 2020
----------------------------------------------------------------


----------------------------------------------------------------
Wiki for the analysis: 
https://wiki.icecube.wisc.edu/index.php/Search_for_neutrinos_from_precursors_and_afterglows_of_GRBs

software used for the analysis: 
csky (revision 176406).

cvmfs environment: 
`/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/setup.sh`
----------------------------------------------------------------


----------------------------------------------------------------
Codes in this repository:

1. Final_analysis.ipynb : The main analysis code. I run this as a jupyter notebook to perform the analysis.
2. BP_distribution.py : The code required to produce the Best Binomial distribution file.  The default produces the result for 1000 scrambled trials.
3. Background_10k_trials_Afterglow.py: The code required to produce prompt+afterglow searches on Background only data for 10,000 trials for each GRB.
4. Background_10k_trials_Precursor.py: The code required to produce Precursor searches on Background only data for 10,000 trials for each GRB.
----------------------------------------------------------------


----------------------------------------------------------------
The main data files required for the analysis are:

1. GRB source list: can be obtained from GRBWeb, simply remove the top heading row, it has been loaded in the main analysis file as 'GRBs_withoutheader.txt'.
2. Afterglow background trials: Can be generated using the script 'Background_10k_trials_Afterglow.py'. This is used to get the p-values from the TS for a 'prompt+afterglow' search.
3. Precursor background trials: Can be generated using the script 'Background_10k_trials_Precursor.py'. This is used to get the p-values from the TS for a 'precursor' search.
4. Distribution of the Best Binomial probabilities (BP_scrambled_bckg.p): This can be generated using the script 'BP_distribution.py'. This is used to obtain the final post-trial pvalue.
----------------------------------------------------------------

----------------------------------------------------------------
The main script which executes the complete analysis has been commented at in appropriate points.
In order to repeat the entire analysis, I do the following:

1. Create the py3-v4 environment using cvmfs.
2. Run 'Background_10k_trials_Afterglow.py' and 'Background_10k_trials_Precursor.py'.
3. Run 'BP_distribution.py'.
4. Run the jupyter notebook 'Final_analysis.ipynb'.
----------------------------------------------------------------

----------------------------------------------------------------
Notes:

1. Make sure that the appropriate data files have the correct paths. 
----------------------------------------------------------------


