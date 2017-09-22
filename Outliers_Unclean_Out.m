
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: Making a subset of UNCLEAN data without outliers (outlier extraction)
% to see outlier criteria see 'Outliers_Clean_Out.m'


%%%%%%%%% Info %%%%%%%%%%%
%%% connectivity - se README
%%% Unclean1 - Vector to subset unclean data relevant to the study

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bring the relevant subsets in
addpath /UndergraduateThesis
Matrices_Out
Outliers_Clean_Out


%%%%%%%%%%%%%%%% Creating variables for indentifiying outliers  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create matrix with relevan Unclean data only
conn_unclean = connectivity(:,:,:, Unclean1);

%%%% Create matrix 'S' storing mean FA (1), mean NOS(2), mean prevalence of
% connections(3) and of disconnections(4)

%--- Storing average NOS, FA per subject
m = size(conn_unclean,4);
S = zeros(m,2);
for i = 1:m
    S(i,1) = mean(nonzeros(conn_unclean(:,:,1,i))); %storing mean NOS
    S(i,2) = mean(nonzeros(conn_unclean(:,:,3,i))); %storing mean FA
end

P = mean(conn_unclean(:,:,1,:)>0, 4);
for i = 1:m %Storing mean P of connections
  B= conn_unclean(:,:,1,i)>0 ;
   S(i,3) = mean(P(B==1));
end

for i = 1:m %Mean prevalence of disconnection
    B = conn_unclean(:,:,1,i)>0;
    C = B + eye(size(B));
    S(i,4) = mean(P(C==0));
end

%%%%%%%%%%% Dividing into quartiles - Obtaining outliers %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k= [25,50,75]; %percentiles wanted
percentiles = prctile(S,k);
IQR= iqr(S); %IQR of the values in columns of S


%--- Collect outlier thresholds for each of the stored values in S
for i=1:4
    outlier_threshold(i)= IQR(i)*1.5;
end

%%%%%%%%%%%%%%%% UPPER Threshold Outliers %%%%%%%%%%%%%%%%%%%%

%---UPPER Threshold (Q3+Threshold)---
for i= 1:4
    Q3s=percentiles(3,i);
    UpperLim(i)= Q3s + outlier_threshold(i);
end

%--- Observing if there are outliers
S(S(:,1)>UpperLim(1),1) %NOS
S(S(:,2)>UpperLim(2),2) %FA
S(S(:,3)>UpperLim(3),3) %Conn Prevalence
S(S(:,4)>UpperLim(4),4) %Disconn Prevalence


%--- Linking the outlier values to the ID in 'subjects'
IndexNumberUnclean = zeros(size(Unclean1,2), 1);
for i = 1:size(Unclean1,2)
    IndexNumberUnclean(i)= Unclean1(i);
end

for i = 1:size(Unclean1,2)
    S(i,5)= IndexNumberUnclean(i);
end

%--- Obtaining the indexes of outliers
S(S(:,2)>UpperLim(2), 5)
S(S(:,4)>UpperLim(4), 5)

UpperOutUn = [284; 792;795]
Unclean2 = Unclean1(~ismember(Unclean1, UpperOutUn)); %subset goal


%%%%%%%%%%%%%%%% LOWER Threshold Outliers %%%%%%%%%%%%%%%%%%%%

%---Lower Threshold (Q1-Threshold)---
%Select first quartile
for i= 1:4
    Q1s=percentiles(1,i);
    LowerLim(i)= Q1s - outlier_threshold(i);
end

%for all columns of S
S(S(:,1)<LowerLim(1),1)
S(S(:,2)<LowerLim(2),2)%Problem
S(S(:,3)<LowerLim(3),3)
S(S(:,4)<LowerLim(4),4)

%Index of the outlier
S(S(:,2)<LowerLim(2),5)
LowerOutUn=[118];

Unclean3 = Unclean2(~ismember(Unclean2, LowerOutUn)); %subset goal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% RECRIPROCAL TAKE OUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PURPOSE: Taking the outliers of Unclean from Clean and of Clean from Unclean.

%%%%%%%%%%%%%%%%%%% Removing Unclean outliers FROM CLEAN %%%%%%%%%%%%%%%%%%%%%%

%--- Find the corresponding indices in Clean
OutUn = [UpperOutUn; LowerOutUn];

for i= 1:size(Clean2,2)
    found = 0;
    Cl_name= subjects{Clean2(i)}(1:12);

    for y= 1:size(OutUn,1)
        Uncl_Out_name = subjects{OutUn(y)}(1:12);
    if strcmp(Cl_name, Uncl_Out_name)== 1
             found = 1;
             x = y
             break;
      end
     end
     if found == 1
         display(Cl_name);
         display(Clean2(i));
         CleanIndexOut(y)= Clean2(i)
     end
end

%FINAL Subset CLEAN
Clean3 = Clean2(~ismember(Clean2, CleanIndexOut));

%%%%%%%%%%%%%%%%% Removing clean outliers FROM UNCLEAN %%%%%%%%%%%%%%%%%%%%%%

OutCl = [UpperOut; LowerOut];

for i= 1:size(Unclean3,2)
    found = 0;
    Un_name= subjects{Unclean3(i)}(1:12);

    for y= 1:size(OutCl,1)
        Clean_Out_name = subjects{OutCl(y)}(1:12);
    if strcmp(Un_name, Clean_Out_name)== 1
             found = 1;
             %x = y
             break;
      end
     end
     if found == 1
         display(Un_name);
         display(Unclean3(i));
        UncleanIndexOut(y)= Unclean3(i) % I think it makes the 1st one zero in here
        %Because it does not find a match wihtin Unclean (it has been
        %removed already)
     end
end

%Note: It suggests only two names to take out, this is because the first one was
%also an outlier within Unclean, so it was already removed in here.

Unclean4 = Unclean3(~ismember(Unclean3, UncleanIndexOut));

%--- Check that subjects are odered in corresponding order within the two FINAL
% subsets
Name_Check= [subjects(Unclean4)', subjects(Clean3)']
%CONFIRMED!
