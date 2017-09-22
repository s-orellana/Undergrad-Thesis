
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:
% 1) Create two vectors allowing to subset clean and unclean data relevant
% to the study
% 2) Make a sub-subset of Unclean that only includes subjects who have both
% clean and unclean data.

%%%%%%%%% Info %%%%%%%%%%%
%%%% subjects - cell array containing ID number of each subject.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Import data
addpath /UndergraduateThesis
load connectivity_dti_aparc.mat

%---Onbtain indexes of the relevant clean and unclea data
Clean = find(~cellfun(@isempty, strfind(subjects, 'm1_c')));
Unclean= find(~cellfun(@isempty, strfind(subjects, 'm1_uc')));


%--- Find IDs of subjects without clean data
 for x = 1:size(Unclean,2)
     found = 0;
     Unclean_name = subjects{Unclean(x)}(1:12);
     for y = 1:size(Clean,2)
         Clean_name = subjects{Clean(y)}(1:12);
         if strcmp(Unclean_name,Clean_name)==1
             found = 1;
             break;
         end
     end
     if found == 0
         display(Unclean_name);
     end
 end

%--- Obtain indices of subjects without clean data
 out1= find(~cellfun(@isempty, strfind(subjects, '16691cu3112d_s_m1_uc')));
 out2= find(~cellfun(@isempty, strfind(subjects, '251cd7u46d39_s_m1_uc')));
 out3= find(~cellfun(@isempty, strfind(subjects, '4b09d6u24c4d_s_m1_uc')));
 out4= find(~cellfun(@isempty, strfind(subjects, '8076adu47184_s_m1_uc')));
 out = [out1, out2, out3, out4];



 %---Create a subset of "Unclean" that does not contain the subject above
 Unclean1 = Unclean(~ismember(Unclean, out));
