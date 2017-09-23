
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT: Mean diffusivity - MD (6).
% PURPOSE: Compute NORMALIZED local Clustering Coefficients for each network
% and compare between clean and unclean.


%%%%%%%%% Info %%%%%%%%%%%
%%% MD_Variables - Output from Global Analyses
%%% C_local_rand_cl/un -  82x1000x275 Matrix.  DIM1: Nodes, DIM2: the local
%%%                       clustering coefficient of 'jth' random network for
%%%                       the subject specified by DIM3.
%%% MD_Node_Out_A - Output from local NON-normalized analyses
%%% LocalClustering_cl/un - 82x275 Matrix. Each row stores the local clustering
%%%                         coefficient of the 'jth' subject (DIM2).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%--- Add paths
addpath /UndergraduateThesis
addpath /UndergraduateThesis/RUN
addpath /UndergraduateThesis/RUN/'Brain Connectivity Toolbox'/
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/Output_Global
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/Output_Local
addpath /UndergraduateThesis/RUN/2_Node_Level_Differences/Normalized_Local_Script/

%--- Introduce connectivity Matrices
Outliers_Unclean_Out
clc
clearvars -except 'Clean3' 'Unclean4' 'connectivity' 'regionDescriptions' 'ROIs' 'subjects' 'weightDescriptions'
conn_clean= connectivity(:,:,:, Clean3);
conn_unclean= connectivity(:,:,:, Unclean4);
m = size(conn_clean(:,:,6,:),4);
n = size(conn_unclean(:,:,6,:),4);


%--- Load data generated in GLOBAL analyses
load ('MD_Variables',  'C_local_rand_cl', 'C_local_rand_un')

%--- Loading data generated in LOCAL non-normalized anayses
load ('MD_Node_Out_A', 'LocalClustering_cl', 'LocalClustering_un')

%%%%%%%%%%%  Normalized Local Clustering Coefficient %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--- Rename (for simplicity) matrix with the local clustering coefficients for
% the random networks of each subject. - Call it DOG
DOG_cl= C_local_rand_cl;
DOG_un= C_local_rand_un;

%--- Result from below: In each entry of the third DIM, I now have a single
% vector (instead of a whole matrix per subject), with the average local
% clustering coefficient for all the random networks of the kth subject.
DOG2_cl = mean(DOG_cl,2);
DOG2_un = mean(DOG_un,2);


%--- Take each vector from the 3rd DIM and store each subject next to
% each other in a matrix.
LCCRand_CL = zeros(82, 275);
for i= 1:size(DOG2_cl,3)
  LCCRand_CL(:,i)= DOG2_cl(:,1,i);%Local Clust Coeff Random Networks - For Clean
end


LCCRand_UN = zeros(82, 275);
for i= 1:size(DOG2_un,3)
  LCCRand_UN(:,i)= DOG2_un(:,1,i);%Local Clust Coeff Random Networks - For Clean
end


%--- Normalized local clusterings
LocalC_Normalized_cl = LocalClustering_cl./LCCRand_CL;
LocalC_Normalized_un = LocalClustering_un./LCCRand_UN;


%----CHECK UP
% I performed this check up to see where the problem with the NaN's was
%[LocalC_Normalized_cl(79,:)', LocalClustering_cl(79,:)', LCCRand_CL(79,:)']
%[LocalC_Normalized_cl(46,:)', LocalClustering_cl(46,:)', LCCRand_CL(46,:)']
'

%making NaNs = 0
LocalC_Normalized_cl(isnan(LocalC_Normalized_cl))=0;
LocalC_Normalized_un(isnan(LocalC_Normalized_un))=0;

%--- ttests
LCCNorm_Comparison = zeros(82, 9);
for i = 1: 82
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(LocalC_Normalized_cl(i,:), LocalC_Normalized_un(i,:), 'alpha', 0.00060975609);
LCCNorm_Comparison(i,1) =  H;
LCCNorm_Comparison(i,2) =  P;
LCCNorm_Comparison(i,3) = STATS.tstat;
LCCNorm_Comparison(i,4) = mean(LocalC_Normalized_cl(i,:));
LCCNorm_Comparison(i,5) = std(LocalC_Normalized_cl(i,:));
LCCNorm_Comparison(i,6) = mean(LocalC_Normalized_un(i,:));
LCCNorm_Comparison(i,7) = std(LocalC_Normalized_un(i,:));
LCCNorm_Comparison(i,8) = STATS.sd;
LCCNorm_Comparison(i,9) = i; %gives you the number of the node
end


%%%%%%%%% WITHOUT BONFERRONI

LCCNorm_Comparison_bon = zeros(82, 9);
for i = 1: 82
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(LocalC_Normalized_cl(i,:), LocalC_Normalized_un(i,:));
LCCNorm_Comparison_bon(i,1) =  H;
LCCNorm_Comparison_bon(i,2) =  P;
LCCNorm_Comparison_bon(i,3) = STATS.tstat;
LCCNorm_Comparison_bon(i,4) = mean(LocalC_Normalized_cl(i,:));
LCCNorm_Comparison_bon(i,5) = std(LocalC_Normalized_cl(i,:));
LCCNorm_Comparison_bon(i,6) = mean(LocalC_Normalized_un(i,:));
LCCNorm_Comparison_bon(i,7) = std(LocalC_Normalized_un(i,:));
LCCNorm_Comparison_bon(i,8) = STATS.sd;
LCCNorm_Comparison_bon(i,9) = i; %gives you the number of the node
end

save('MD_Clust_Norm_Local_Output')
