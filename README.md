# Undergraduate Thesis

This repository contains the code implemented in some of the data analysis steps of my bachelor's thesis at University College Utrecht, titled: '*The Influence of FreeSurfer’s Manual Editing Options on Structural Brain Network Reconstruction*'. 
Written in Spring 2017 at the [UMC Utrecht](http://www.umcutrecht.nl/nl/).

## Introduction
In the resconstruction of structural brain networks, a common choice of software for cortical parcellation and subsequent node defintion is Freesurfer's (FS) processing pipeline (http://surfer.nmr.mgh.harvard.edu/). Despite its high reliability, mistakes can still be made by FS, and therefore the introduction of manual edits is allowed. 

The goal of my thesis was to determine whether the topology of structural brain networks which were reconstructed with a fully automated FS pipeline would differ significantly from those reconstructed in a semi-automated fashion (i.e. the introduction of manual edits was allowed). The code written for a section of the data analysis is presented in this repository.


## Data Analysis
The data for this code is provided by the file 'connectivity_dti_aparc.mat'. Amongst other things, it contains connectivity matrices (describing structural brain networks) generated with both edited and unedited Freesurfer data.


### Structure of the input data
The file connectivity_dti_aparc.mat, unavailable here, contains the following:

1) A four dimensional (82x82x8x900) matrix called 'connectivity' organized as follows:
    - DIMs 1 and 2 pertain to the nodes of the network.
    - DIM 3 gives the weight modality for the connections (i.e. value of matrix entries). There are seven possible modalities. I use the following: (1) NOS, (2) Fiber distance, (3) FA, (6) mean diffusivity, (7) streamline volume density.
    - DIM 4 indexes subjects. Edited and Unedited data of the same subject are stored as if they were different subjects.

2) A cell array titled 'subjects' containing the IDs of each subject. Their order corresponds to entries in DIM4 of the matrix. The IDs also codify whether the data in that particular entry is clean or unclean.

(Note: the only contents described are those relevant to the code)


### General steps

**Step 1** : Extracts the data that is eligible for the study (i.e. making sure only subjects with both edited and unedited data are included). See 'Matrices_Out' 'Outliers_Clean_Out' ' Outliers_Unclean_Out'

**Step 2**: Computes *global* graph theory measures and compares them in edited vs. unedited networks. See scripts in 'RUN/1_Global_Comparison'.

**Step 3**: Computes *local* graph theory measures and compares them in edited vs. unedited networks. For non-normalized metrics see scripts in 'RUN/2_Node_Level_Differences'. For normalized metrics see scripts in 'RUN/2_Node_Level_Differences/ Normalized_Local_Script'

Note: Analyses were done for both weighted and unweighted networks. This repository only contains code for the weighted analyses. 


## Notes
Some sections of the code in this repository were written during the '10K in a Day' workshop by the [Dutch Connectome Lab](http://www.dutchconnectomelab.nl/)

- **Clean** refers to data that suffered *manual editing* in FreeSurfer. **Unclean** refers to data that received *no manual editing* in the FS pipeline.
- Scripts usually call for the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) (Rubinov & Sporns, 2010).

## References
Rubinov M, Sporns O (2010) NeuroImage 52:1059-69
