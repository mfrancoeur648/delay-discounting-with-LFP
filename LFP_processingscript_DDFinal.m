%% initial code to extract and process data from Open Ephys (XDF) files (behavior, LFP). 
% 
% DDFinal modification 

%
clear all;

xdfdir = 'C:\Users\Neatlabs Data Engine\Documents\MATLAB\load_xdffunc';
xdfdir2 = 'C:\Users\Neatlabs Data Engine\Documents\MATLAB\xdf-Matlab-master';
addpath(xdfdir)
addpath(xdfdir2)

rdir='E:\DelayDiscounting_Final_LFP_raw\Pharmacology Study\';  % point to directory where raw data is kept. 

%Assuming current organization of raw data is based on data acquiring. Will
%resort/save data for each animal based on data acquired.

%% switched paradigm from 5 Sec LR delay to 500ms LR delay after..(

ddir=dir(rdir);            % variable for the list of folders that need to be processed.
ddir(1:2)=[];
rerunany=1; % if you want to rerun any part of the analysis default 1
rerun.beh=1; % if you want to re-run data extraction; otherwise will skip files where extracted default 1
rerun.ERP=0; % if you want to re-run ERP script; otherwise will skip files where extracted
rerun.ERSP=0; % if you want to re-run ERSP script(long time); otherwise will skip files where extracted
cd(rdir);
%par=load ('preprocparams.mat');
par=0;
%%

drugType = ["C","M","S","BO"]; % manually designate to store data separately for each drug type

tic
for d = 1:length(drugType)
    sdir=char(strcat('E:\dd_final_processed\DrugTypeSorted\',drugType(d),'\'));    %point to directory where processed data will be kept.
    if ~exist(sdir,'dir')
        mkdir(sdir)
    end
    for i =1:length(ddir)  % iterate through all folders in raw data dir
        disp(ddir(i).name);
        if(startsWith(ddir(i).name,['DD_08312022'])) % edit to run for specfiic folder 
            cd([rdir '\' ddir(i).name]);
            adir=dir;
            adir(1:2)=[];
            for j = 1:length(adir) % iterate through all animals recorded on that date
                an=adir(j).name;
                ani = extractBefore(an,'_'); % store animal name
                if(contains(an,drugType(d)))
                    if(~contains(an,'B')) % check to see if B program was run
                        savedir=[sdir ani '\' ddir(i).name '\'];
                        anf=[adir(j).folder '\' an];
                        if (~exist('anf','dir') ||  rerunany==1)
                            extractLFPbeh_DD2to10(anf,savedir,par,rerun);
                            disp(an);
                            disp('entered');
                        end  
                    elseif(contains(an,'B')) % run on B program if included in name
                        savedir=[sdir ani 'B\' ddir(i).name '\'];
                        anf=[adir(j).folder '\' an];
                        if (~exist('anf','dir') ||  rerunany==1)
                            extractLFPbeh_DD2to10B(anf,savedir,par,rerun);
                            disp(an);
                            disp('entered');
                            disp('B file found')
                        end  
                    end
                end
            end 
            disp(i);
        end
    end
end
toc
disp('Raw data xdf processing complete');

