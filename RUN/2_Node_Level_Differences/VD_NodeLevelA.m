
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT: Streamline Volume Density - VD (7)
% PURPOSE: Compute measures on the local topology of the networks and
% compare them. NON-normalized!


%%%%%%%%% Info %%%%%%%%%%%
%%% NOS_Node_Out_A - stores the outputs of this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%--- Inport data, add paths
addpath /UndergraduateThesis
addpath /UndergraduateThesis/2_Node_Level_Differences
Outliers_Unclean_Out
clc
addpath /UndergraduateThesis/'Brain Connectivity Toolbox'/
clearvars -except 'Clean3' 'Unclean4' 'connectivity' 'regionDescriptions' 'ROIs' 'subjects' 'weightDescriptions'

%--- creating connectivity matrices
conn_clean= connectivity(:,:,:, Clean3);
conn_unclean= connectivity(:,:,:, Unclean4);

m = size(conn_clean(:,:,7,:),4);
n = size(conn_unclean(:,:,7,:),4);

%%%%%%%%%%%%%%%%%%%%%%%% NODE STRENGTH %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%----------Clean-------------

NSM_CL= zeros(size(conn_clean,1), size(conn_clean,4)); %Node Strengths Matrix
%MAT where each jth column represents a subject, each ith row represents
%the strength of the node.

for i = 1:m
A = conn_clean(:,:,7,i);
NSM_CL(:,i)= sum(A,2);
end


%----------Unclean-------------
NSM_UN= zeros(size(conn_unclean,1), size(conn_unclean,4));

for i = 1:m
A = conn_unclean(:,:,7,i);
NSM_UN(:,i)= sum(A,2);
end


%--- ttests
NSM_Comparison = zeros(size(NSM_UN,1), 9);
for i = 1: size(NSM_UN,1)
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(NSM_CL(i,:), NSM_UN(i,:), 'alpha', 0.00060975609);
 NSM_Comparison(i,1) =  H;
 NSM_Comparison(i,2) =  P;
 NSM_Comparison(i,3) = STATS.tstat;
 NSM_Comparison(i,4) = mean(NSM_CL(i,:));
 NSM_Comparison(i,5) = std(NSM_CL(i,:));
 NSM_Comparison(i,6) = mean(NSM_UN(i,:));
 NSM_Comparison(i,7) = std(NSM_UN(i,:));
 NSM_Comparison(i,8) = STATS.sd;
 NSM_Comparison(i,9)= i; %number of the node
end

%{
%What is stored in each column of NSM_Comparison
PER ROW: the row represents the ith node. Each column is (usually) the result
of the comparison between the strengths of the ith node in clean v unclean
(1) = H value of the t-test comparison
(2) = P value of the t-test
(3) = t-statistic of the t-test
(4) = mean of the strengths for the ith node, in CLEAN
(5) = STD of the strengths, in CLEAN
(6) = mean of the strengths,  in UNCLEAN
(7) = STD of the strengths, in UNCLEAN
(8) = STD of the difference (?)
(9) = number of the node
%}


%%%%%%%%%%%%%%%%%%LOCAL CLUSTERING COEFFICIENT %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Non-Normalized %%%%%%%%%%%

%----Clean
LocalClustering_cl=(zeros(size(conn_clean,1), size(conn_clean,4)));
for i = 1:m
A = conn_clean(:,:,7,i);
LocalClustering_cl(:,i) = clustering_coef_wu(A);
end

%----Unclean
LocalClustering_un=(zeros(size(conn_unclean,1), size(conn_unclean,4)));
for i = 1:m
A = conn_unclean(:,:,7,i);
LocalClustering_un(:,i) = clustering_coef_wu(A);
end

%--- ttests
LClus_Comparison = zeros(82, 9);
for i = 1: 82
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(LocalClustering_cl(i,:), LocalClustering_un(i,:), 'alpha', 0.00060975609);

LClus_Comparison(i,1) =  H;
LClus_Comparison(i,2) =  P;
LClus_Comparison(i,3) = STATS.tstat;
LClus_Comparison(i,4) = mean(LocalClustering_cl(i,:));
LClus_Comparison(i,5) = std(LocalClustering_cl(i,:));
LClus_Comparison(i,6) = mean(LocalClustering_un(i,:));
LClus_Comparison(i,7) = std(LocalClustering_un(i,:));
LClus_Comparison(i,8) = STATS.sd;
LClus_Comparison(i,9) = i; %gives you the number of the node
end

%{
%What is stored in each column of NSM_Comparison
PER ROW: the row represents the ith node. Each column is (usually) the result
of the comparison between the clustering coeffs of the ith node in clean v unclean
(1) = H value of the t-test comparison
(2) = P value of the t-test
(3) = t-statistic of the t-test
(4) = mean of the clustering coeffs for the ith node, in CLEAN
(5) = STD of the clustering coeffs, in CLEAN
(6) = mean of the clustering coeffs,  in UNCLEAN
(7) = STD of the clustering coeffs, in UNCLEAN
(8) = STD of the difference (?)
(9) = number of the node
%}




%%%%%%%%%%%%%%%%%%%%%%%% LOCAL EFFICIENCY %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Non-normalized %%%%%%%%

%----Clean
LocalEfficiency_cl=(zeros(size(conn_clean,1), size(conn_clean,4)));
for i = 1:m
A = conn_clean(:,:,7,i);
LocalEfficiency_cl(:,i) = efficiency_wei(A,1);
WeAreAt_LocalE_CL =  i
WHERE = 'VD'
end


%----Unclean
LocalEfficiency_un=(zeros(size(conn_unclean,1), size(conn_unclean,4)));
for i = 1:m
A = conn_unclean(:,:,7,i);
LocalEfficiency_un(:,i) = efficiency_wei(A,1);
WeAreAt_LocalE_UN =  i
WHERE = 'VD'
end


%--- ttests
LEfficiency_Comparison = zeros(82, 9);
for i = 1: 82
    %bonferroni correction - P= 0.05/82
[H,P,CI,STATS]= ttest(LocalEfficiency_cl(i,:), LocalEfficiency_un(i,:), 'alpha', 0.00060975609);
LEfficiency_Comparison(i,1) =  H;
LEfficiency_Comparison(i,2) =  P;
LEfficiency_Comparison(i,3) = STATS.tstat;
LEfficiency_Comparison(i,4) = mean(LocalEfficiency_cl(i,:));
LEfficiency_Comparison(i,5) = std(LocalEfficiency_cl(i,:));
LEfficiency_Comparison(i,6) = mean(LocalEfficiency_un(i,:));
LEfficiency_Comparison(i,7) = std(LocalEfficiency_un(i,:));
LEfficiency_Comparison(i,8) = STATS.sd;
LEfficiency_Comparison(i,9) = i; %gives you the number of the node
end


save (UndergraduateThesis/2_Node_Level_Differences/Output_Local/'VD_Node_Out_A.mat')
