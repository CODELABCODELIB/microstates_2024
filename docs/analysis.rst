Methods
=======
This code identifies microstates in EEG data while participants are engaged in smartphone behavior. The relationship between the microstates and behavior is explored.

Pre-processing EEG
------------------
1. Loaded the EEG struct and included participants that met the following criterias:

 - Had EEG data
 - Contains smartphone data
 - Has been Aligned to Smartphone data
 - If the participant was a curfew participant, we only included the first measurement file
 - Attys is false
 - Measurement did not contain a saving error

2. Cleaned the EEG data by running gettechnicallcleanEEG

 - Removed blinks according to ICA
 - Interpolated missing channels
 - Re-referenced data to the average channel
 - Highpass filter up to 45 Hz
 - Lowpass filter from 0.5 Hz

Identify prototypical microstates per participant with toolbox[1].
==================================================================
1.  Select EEG of participants while they are using their smartphone and identify GFP peaks. Datatype is spontaneous as the data is not timelocked. Minimum Distance between GFP peaks is 10 ms and the clustering will be based on 1000 randomly selected GFP peaks.

2.  No norrmalization of the data is performed (?)

3. Cluster the GFP peaks with modified k-means algorithm. Identifying 2 to 10 microstates. Repeating the algorithm with 10 repetitions, and capping the number of iterations at 1000. Convergence threshold is set to 1e-6. Note that modified k-means is polarity invariant.

4. Select the optimal number of clusters (between 2 to 10) by optimal Global explained variance (GEV) 

5. Fit the protypical maps identified in step 4 to the overal data (every timepoint)

6. No smoothing is performed (?)

7. Calculate statistics per participants:

 - TODO

[1] https://archive.compute.dtu.dk/files/public/users/atpo/Microstate. Poulsen, A. T., Pedroni, A., Langer, N., & Hansen, L. K. (2018). Microstate EEGlab toolbox: An introductory guide. BioRxiv, 289850.

Joint probability interval distribution (JID-state)
===================================================
1. Assign each touchscreen interaction to nearest 'bin/pixel'. 
2. Select the active microstate at the time of the interaction
3. Repeat for each interaction
4. Create a JID-state for each microstate per participant

Population level effects 
========================
Cluster participants protopical maps with k-means to find similar clusters across participants 

1. Perform cross-validation to identify optimum number of 'k'

2. TODO
