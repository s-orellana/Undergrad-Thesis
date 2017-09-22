
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:
% Making a subset of CLEAN data without outliers (outlier extraction)
% NOTE: Outliers are defined to be subject whose average FA or NOS values are
% above 'Q3 + 1.5* IQR' or below 'Q1 + 1.5* IQR'; IQR= interquatile range
% Q3= boundary between 3rd and 4ht quartiles; Q1= boundary between 1st and 2nd
% quartile.

%%%%%%%%% Info %%%%%%%%%%%
%%%% Clean - vector with indexes of clean data
%%%% connectivity - see README

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bring the relevant subset in
%addpath /UndergraduateThesis
Matrices_Out

%%%%%%%%%%%%%%%% Creating variables for indentifiying outliers  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create matrix with Clean data only
conn_clean = connectivity(:,:,:, Clean);


%%%% Create matrix 'S' storing mean FA (1), mean NOS(2), mean prevalence of
% connections(3) and of disconnections(4)

%--- Storing average NOS, FA per subject
m = size(conn_clean,4);
S = zeros(m,2);
for i = 1:m
    S(i,1) = mean(nonzeros(conn_clean(:,:,1,i))); %storing mean NOS
    S(i,2) = mean(nonzeros(conn_clean(:,:,3,i))); %storing mean FA
end

P = mean(conn_clean(:,:,1,:)>0, 4);
for i = 1:m %Storing mean P of connections
    B = conn_clean(:,:,1,i)>0 ;
   S(i,3) = mean(P(B==1));
end

for i = 1:m %Mean prevalence of disconnection
    B = conn_clean(:,:,1,i)>0;
    C = B + eye(size(B));
    S(i,4) = mean(P(C==0));
end

%%%%%%%%%%% Dividing into quartiles - Obtaining outliers %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k= [25,50,75]; %percentiles wanted
percentiles = prctile(S,k); %matrix of percentiles (rows percentiles, columns = columns of S)
IQR= iqr(S); %IQR of the values in columns of S


%--- Collect outlier thresholds for each of the stored values in S
for i=1:4
    outlier_threshold(i)= IQR(i)*1.5;
end

%%%%%%%%%%%%%%%% UPPER Threshold Outliers %%%%%%%%%%%%%%%%%%%%

%---UPPER Threshold (Q3+Threshold)---
for i= 1:4
    Q3s=percentiles(3,i); %Select third quartile
    UpperLim(i)= Q3s + outlier_threshold(i);
end

%--- Observing if there are outliers
%for NOS
S(S(:,1)>UpperLim(1),1);
%for FA
S(S(:,2)>UpperLim(2),2)
% Conn Prevalence
S(S(:,3)>UpperLim(3),3);
% Disconn Prevalence
S(S(:,4)>UpperLim(4),4)

%--- Linking the outlier values to the ID in 'subjects'

%List the indexes of Clean
IndexNumberClean = zeros(size(Clean,2), 1);
for i = 1:size(Clean,2)
    IndexNumberClean(i)= Clean(i);
end

%Adding the subject index (corresponding to its name in 'subjects' & position
% in 'connectivity') to 'S'
for i = 1:size(Clean,2)
    S(i,5)= IndexNumberClean(i);
end

%--- Obtaining the indexes of the ourliers
S(S(:,2)>UpperLim(2), 5)
S(S(:,4)>UpperLim(4), 5)

UpperOut= [117; 678]; %Making a vector of exclusions with the indexes
Clean1 = Clean(~ismember(Clean, UpperOut)); %Subset goal


%%%%%%%%%%%%%%%% LOWER Threshold Outliers %%%%%%%%%%%%%%%%%%%%

%--- LOWER Threshold (Q1-Threshold)---
for i= 1:4 %Select first quartile
    Q1s=percentiles(1,i);
    LowerLim(i)= Q1s - outlier_threshold(i);
end

%for all columns of S
S(S(:,1)<LowerLim(1),1)
S(S(:,2)<LowerLim(2),2)
S(S(:,3)<LowerLim(3),3)
S(S(:,4)<LowerLim(4),4)

%Index of the outlier
S(S(:,4)<LowerLim(4),5)
LowerOut=[730];

%subset
Clean2 = Clean1(~ismember(Clean1, LowerOut));
