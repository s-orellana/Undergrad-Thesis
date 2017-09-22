
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT: Fiber Distance - FD (2).
% PURPOSE:
% 1) Compute measures on the global topology of the networks and
% compare them.
% 2) Compute the normalized measures when relevant
% 3) Store created random networks for each subject

%%%%%%%%% Info %%%%%%%%%%%
%%% Clean3: Vector to subset relevant clean data from 'connectivity'
%%% Unclean4: Vector to subset relevant unclean data from 'connectivity'
%%% FD_Variables: file saving the output of this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%Bring data in
addpath /UndergraduateThesis
Outliers_Unclean_Out
clc
addpath /UndergraduateThesis/1_Global_Comparison
addpath /UndergraduateThesis/'Brain Connectivity Toolbox'/

%--- leave only relevant variables in
clearvars -except 'Clean3' 'Unclean4' 'connectivity' 'regionDescriptions' 'ROIs' 'subjects' 'weightDescriptions'

%--- make separate matrices
conn_clean= connectivity(:,:,:, Clean3);
conn_unclean= connectivity(:,:,:, Unclean4);

m = size(conn_clean(:,:,2,:),4);
n = size(conn_unclean(:,:,2,:),4);

%-------------TOTAL WEIGHT (sum of all weights) per subject network (1)---------
%-----------------------------------------------------------------------------
totalW_cl = zeros(m,1);
totalW_un = zeros(n,1);

for i = 1:m
    totalW_cl(i)=(sum(sum(conn_clean(:,:,2,i))));
end

for i = 1:m
    totalW_un(i) = (sum(sum(conn_unclean(:,:,2,i))));
end

%
%Test
Mean_totalW_cl = mean(totalW_cl);
Std_totalW_cl = std(totalW_cl);
Mean_totalW_un = mean(totalW_un);
Std_totalW_un = std(totalW_un);
[h1,p1,ci1,stats1] = ttest(totalW_cl, totalW_un);%not stored as a variable


%------------- Average Connection Weight  per subject network  (2)----------------
%---------------------------------------------------------------------------------
AvConnW_cl= zeros(m,1);
AvConnW_un= zeros(n,1);

for i = 1:m
   AvConnW_cl(i)=  mean(nonzeros(conn_clean(:,:,2,i))); %squareform
end

for i = 1:n
    AvConnW_un(i) = mean(nonzeros(conn_unclean(:,:,2,i)));
end

Mean_AvConnW_cl = mean(AvConnW_cl);
Std_AvConnW_cl = std(AvConnW_cl);
Mean_AvConnW_un = mean(AvConnW_un);
Std_AvConnW_un = std(AvConnW_un);

[h2,p2,ci2,stats2] = ttest(AvConnW_cl, AvConnW_un);

%--------------------------Clustering Coefficient ------------------------
%------------------------------------------------------------------------

%-----------------------------Clean--------------------------
C_CL = zeros(m,1);
C_CL_rand_pre= zeros(m,1000);
C_local_rand_cl= zeros(82,1000,m); %DIM2= number of random networks

for i= 1:m
    A = conn_clean(:,:,2,i);
    C_CL(i) = mean(clustering_coef_wu(A));
    for j= 1:1000
        B =randmio_und(A,10);%obtain a single random network
        C_CL_rand_pre(i,j) = mean(clustering_coef_wu(B));%AA
        C_local_rand_cl(:,j,i)= clustering_coef_wu(B); %BB
    end
end


%AA: For the ith subject(and row) we store in the "jth" column the global
%clustering coefficient of the jth random network -being computed for that
%"ith" subject.

%BB: 3D Matrix - The third DIM is for a subject, rows are the nodes, each
%column is filled with the LOCAL clustering coefficients of the "Jth"
%random network of the "Kth" subject.

C_CL_rand= mean(C_CL_rand_pre,2); %C random for clean
C_CL_normalized = C_CL ./ C_CL_rand;


%------------------------------Unclean--------------------------
C_UN = zeros(n,1);
C_UN_rand_pre= zeros(n,1000);
C_local_rand_un= zeros(82,1000,n); %DIM2= number of random networks

for i= 1:n
    A = conn_unclean(:,:,2,i);
    C_UN(i) = mean(clustering_coef_wu(A));
    for j= 1:1000
        B =randmio_und(A,10);%obtain a single random network
        C_UN_rand_pre(i,j) = mean(clustering_coef_wu(B));
        C_local_rand_un(:,j,i)= clustering_coef_wu(B);
    end
end


C_UN_rand= mean(C_UN_rand_pre,2); %C random for clean  %CHECK HERE
C_UN_normalized = C_UN./ C_UN_rand;



%-------Comparison - NON-normalized--------
MeanC_nonnorm_cl = mean(C_CL);
StdC_nonnorm_cl = std(C_CL);
MeanC_nonnorm_un = mean(C_UN);
StdC_nonnorm_un = std(C_UN);

[h3,p3,ci3,stats3] = ttest(C_UN, C_CL);

%-------Comparison -NORMALIZED--------
MeanC_norm_cl = mean(C_CL_normalized);
StdC_norm_cl = std(C_CL_normalized);
MeanC_norm_un = mean(C_UN_normalized);
StdC_norm_un = std(C_UN_normalized);

[h4,p4,ci4,stats4] = ttest(C_UN_normalized, C_CL_normalized);



%---------------------- Efficiency ----------------------------
%--------------------------------------------------------------


%------------------------------Clean--------------------------
E_CL= zeros(m,1);
E_CL_rand_pre = zeros(m,1000);
E_local_rand_cl = zeros (82, 1000, m); %change the 5 to the number of iterations

for i= 1:m
    A = conn_clean(:,:,2,i);
    E_CL(i) = efficiency_wei(A); %Put here the efficiency here
    for j= 1:1000
        B =randmio_und(A,10);
        E_CL_rand_pre(i,j) = efficiency_wei(B);
        E_local_rand_cl(:,j,i)= efficiency_wei(B,1);
    end
end

E_CL_rand= mean(E_CL_rand_pre,2);
E_CL_normalized = E_CL./E_CL_rand; %normalized


%------------------------------Unclean--------------------------

E_UN= zeros(n,1);
E_UN_rand_pre = zeros(n,1000);
E_local_rand_un = zeros (82, 1000, n); %change the 5 to the number of iterations

for i= 1:n
    A = conn_unclean(:,:,2,i);
    E_UN(i) = efficiency_wei(A); %Put here the efficiency here
    for j= 1:1000
        B =randmio_und(A,10);
        E_UN_rand_pre(i,j) = efficiency_wei(B);
        E_local_rand_un(:,j,i)= efficiency_wei(B,1);
    end
end

E_UN_rand= mean(E_UN_rand_pre,2);
E_UN_normalized = E_UN./E_UN_rand;


%-------Comparison - NON-normalized--------
MeanE_nonnorm_cl = mean(E_CL);
StdE_nonnorm_cl = std(E_CL);
MeanE_nonnorm_un = mean(E_UN);
StdE_nonnorm_un = std(E_UN);

[h5,p5,ci5,stats5] = ttest(E_UN, E_CL);

%-------Comparison -NORMALIZED--------
MeanE_norm_cl = mean(E_CL_normalized);
StdE_norm_cl = std(E_CL_normalized);
MeanE_norm_un = mean(E_UN_normalized);
StdE_norm_un = std(E_UN_normalized);

[h6,p6,ci6,stats6] = ttest(E_UN_normalized, E_CL_normalized);

%Save output to be used on step 2 of analysis
save (UndergraduateThesis/2_Node_Level_Differences/Output_Global/'FD s_Variables')
