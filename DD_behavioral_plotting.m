% Delay Discounting Behavioral Data Plotting and Data Storage
% Plots %Hr for each delay length

clearvars;
%% DD beh data across all delay lengths 
clear all; 
close all;

rdir='E:\dd_processed\trType_labeled\'; 
esdir ='E:\DD_PhysProcessed\rewL_StatsData\';
if~exist(esdir,'dir')
   mkdir(esdir);
end 

pdir = 'E:\DD_PhysProcessed\rewL plots\Beh_Plots\';
if~exist(pdir,'dir')
    mkdir(pdir);
end
cd(rdir)
ddir = dir;
ddir(1:2)=[];
ancnt=1;
saveVars=1;
saveFigs=1;
plotFigs=1;
c1=0;
sescnt=0;
totDLcnt = 0;
% switch based on which is running
tic
numDelays = length(ddir);
for i = 1:numDelays % Iterate through all delay lengths 
    disp(ddir(i).name)         
    cd([rdir '\' ddir(i).name]);
    adir = dir;
    adir(1:2) = []; 
    delayL(i) = str2num(extractBefore(ddir(i).name,'msH')); %#ok<ST2NM>
    totDLcnt = length(adir);
    ancnt=1;
    sescnt=0;
    for j = 1:totDLcnt %Iterate through all sessions(across all anis)
        an = adir(j).name;
        if(startsWith(an,'R')&&contains(an,'_DD_'))                
            rat = extractBefore(an,'_DD_');
            if ancnt>1
%                 if(saveVars)
%                     if(exist('D','var'))
%                         savedir = strcat(sdir,num2str(delayL(i)),'msDelay\',rname{ancnt-1},'\');
%                         if~exist(savedir,'dir')
%                            mkdir(savedir);
%                         end 
%                         fname = strcat(rname{ancnt-1},'_rewLbinned.mat');
%                         save(fullfile(char(savedir),char(fname)),'D','-v7.3');
%                         clear D;
%                         disp(rname{ancnt-1});
%                     end
%                 end
                if ~strcmp(rat,rname{ancnt-1})
                    rname{ancnt} = rat;
                    ancnt=ancnt+1;
                    sescnt=0;
                end
            else
                rname{ancnt} = rat;
                ancnt=ancnt+1;
            end
            cd([adir(j).folder '\' an]);
            %if(exist('epochdata.mat','file')) % Only include behavioral data with phys            
                sescnt = sescnt+1;
                D.An_phys{sescnt,1}=rat;
                D.An_phys{sescnt,2}=delayL(i);
                if(exist('extracted_data.mat','file'))
                    load('extracted_data.mat');
                    date = extractBetween(an,'_DD_','_');
                    %load('epochdata2.mat');
                    D.behm= behm;
                    numTr = behm.numHighTr + behm.numLowTr;
                    % store for each animal and session in each DL- stats
                    DL(i).highP{ancnt,sescnt+1} = behm.numHighTr/numTr;
                    DL(i).highP{ancnt,1} = rname{ancnt-1};
                    hp(j) =  behm.numHighTr/numTr;
                end
            %end
        end
    end
    % Mean and SEM across all animals and sessions, store for each DL
    hcP(i) = mean(hp,'omitnan');
    hcP_sem(i) = std(hp,'omitnan')/sqrt(totDLcnt);

    if(saveVars) % Save sheet for each delay length
        filename='DD_behavioral_statsdata.xlsx';
        hcpT = cell2table(DL(i).highP);
        numSess = width(hcpT);
        seslbl = [];
        for s = 1:numSess
            seslbl = [seslbl,'Sesssion ' +num2str(s)];
        end
        hcpT.Properties.VariableNames = seslbl;
        dStr = ddir(i).name;
        writetable(hcpT,fullfile(esdir,filename),'FileType','spreadsheet','Sheet',dStr);
        clear hcpT
        % Remove Extra Blank Sheets
        excelFileName=filename;
        excelFilePath = esdir;
        sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
        % Open Excel file.
        objExcel = actxserver('Excel.Application');
        objExcel.Workbooks.Open(fullfile(excelFilePath, excelFileName)); % Full path is necessary!
        % Delete sheets.
        try
              % Throws an error if the sheets do not exist.
              objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
              objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
              objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
        catch
               % Do nothing.
        end
        % Save, close and clean up.
        objExcel.ActiveWorkbook.Save;
        objExcel.ActiveWorkbook.Close;
        objExcel.Quit;
        objExcel.delete;
    end
end
delays = 1:numDelays;
if(plotFigs)
    figure;
    hold on;
    [sortedD, sIdx] = sort(delayL);
    tmp = hcP(sIdx);
    tmp_sem = hcP_sem(sIdx);
    errorbar(tmp,tmp_sem);
    %     shadedErrorBar(sesvec,tbi_tcp,tbi_tcp_sem,...
%         'lineprops',colvec(2),'patchSaturation',0.33); 
    xlabel('Delay Lengths');
    xticklabels(sortedD);
    xticks(1:numel(sortedD));
    %ylim([0,1]);
    ylabel('High Choice Percentage');    
    flabel='High Choice Percentage across Delay Lengths';
    title(flabel); 
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
    %close all;
end
toc
disp('DD behavioral data plotted and stored');
    



