% Data storage for Stats script DD phys
% Store high and low reward trial types 
% Store as sheets corresponding to delays containing matrix of beta activity in 12 rew
% related electrodes across all animals and sessions
% Version 2 : Store for all frequencies
clearvars;
%% Preprocess store mean across sessions across animals for each delay length
rdir= 'E:\DD_PhysProcessed\rewL_processedDataV2_New2\';
pltdir = 'E:\DD_PhysProcessed\rewL plots\';
if ~exist(pltdir,'dir')
    mkdir(pltdir)
end
cd(rdir);
ddir=dir(rdir);
ddir(1:2)=[];

load('E:\DD_PhysProcessed\times_freqs.mat');
% Set Desired Time Window for plotting 

tw1 = -500;
tw2 = 2500;

twa1 = tw1-30;
twa2 = tw2+30;
tt = t(t>twa1 & t<twa2);
ttvec = find(t==min(tt)):find(t==max(tt));

% Early Time Window
etw1 = 0;
etw2 = 1000;

etwa1 = etw1-30;
etwa2 = etw2+30;
tt_e = t(t>etwa1 & t<etwa2);
etvec = find(t==min(tt_e)):find(t==max(tt_e));
ebinidx1 = min(etvec);
ebinidx2 = max(etvec);

% Late Time Window
ltw1 = 1000;
ltw2 = 2000;

ltwa1 = ltw1-30;
ltwa2 = ltw2+30;
tt_l = t(t>ltwa1 & t<ltwa2);
ltvec = find(t==min(tt_l)):find(t==max(tt_l));
lbinidx1 = min(ltvec);
lbinidx2 = max(ltvec);

% Base Line Time Windo
bstw1 = -2500;
bstw2 = -1500;

bstwa1 = bstw1-30;
bstwa2 = bstw2+30;
tt_b = t(t>bstwa1 & t<bstwa2);
bslvec = find(t==min(tt_b)):find(t==max(tt_b));


% Base Line Time Window 2ms delay 
bstw1 = -2500;
bstw2 = -2000;

bstwa1 = bstw1-30;
bstwa2 = bstw2+30;
tt_b = t(t>bstwa1 & t<bstwa2);
bslvec2 = find(t==min(tt_b)):find(t==max(tt_b));



saveVars=1;
saveFigs=1;
plotFigs=0;
c1=0;
sescnt=0;
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
evec = [28,27,26,25,31,32,29,30,6,5,1,20]; % 10 target regions in same order as Go/Nogo figure 
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
twlab1 = tw1:((tw2-tw1+1)/length(ttvec)):tw2;
% switch based on which is running
ld_cnt=1;
md_cnt=1;
hd_cnt=1;                
tic
for b = 1:6
    for i = 1:length(ddir) % Iterate through all delay lengths
        dcnt = 1; % reset count for each delay
        disp(ddir(i).name)
        cd([rdir '\' ddir(i).name]);
        adir = dir;
        adir(1:2) = [];
        ancnt=0;
        Delay_len = str2num(extractBefore(ddir(i).name,'ms'));
        for j = 1:length(adir) %Iterate through all animals
            ancnt = ancnt+1;
            an = adir(j).name;
            rname= an;
            cd([adir(j).folder '\' an]);
            idir = dir;
            idir(1:2)= [];
            if(endsWith(idir.name,'rewLbinned.mat'))
                load(idir.name)
                pdir = strcat(pltdir,ddir(i).name);
                if ~exist(pdir,'dir')
                    mkdir(pdir);
                end
                % Plot DD rewL power time spectograms for every animal session
                % and across all delay lengths
                if ~isfield(D,'rewL_hr_blc')
                    continue;
                end
                E = size(D.rewL_hr_blc,2);
                numSes = size(D.rewL_hr_blc,1);
                numEs = size(D.rewL_hr_blc,2);
                lscnt = 0;
                mscnt = 0;
                hscnt = 0;
                for s = 1:numSes
                    ecnt=4; %first 4 columns are animal and session labels 
                    % delay len and %hr ch
                    for e = evec
                        if e > numEs % skip sessions without all electrodes
                            continue;
                        end
                        ecnt= ecnt+1;
                        nanvec = NaN(1,200);
                        if(length(cell2mat(D.rewL_hr_blc(s,e,b)))~=200)
                            D.rewL_hr_blc{s,e,b} = nanvec;
                        end
                        if(length(cell2mat(D.rewL_lr_blc(s,e,b)))~=200)
                            D.rewL_lr_blc{s,e,b} = nanvec;
                        end
%                         tmp_h = zscore(D.rewL_hr_blc{s,e,b});
%                         tmp_l = zscore(D.rewL_lr_blc{s,e,b});
                        % Subtract by baseline mean, normalize by dividing
                        % by baseline standard deviation
                        if(Delay_len == 2)
                            tmp_h = (D.rewL_hr_blc{s,e,b}-mean(D.rewL_hr_blc{s,e,b}...
                            (bslvec2)))./std(D.rewL_hr_blc{s,e,b}(bslvec2));
                            tmp_l = (D.rewL_lr_blc{s,e,b}-mean(D.rewL_lr_blc{s,e,b}...
                                (bslvec2)))./std(D.rewL_lr_blc{s,e,b}(bslvec2));
                            tmp_df = tmp_h-tmp_l;

                        else
                            tmp_h = (D.rewL_hr_blc{s,e,b}-mean(D.rewL_hr_blc{s,e,b}...
                                (bslvec)))./std(D.rewL_hr_blc{s,e,b}(bslvec));
                            tmp_l = (D.rewL_lr_blc{s,e,b}-mean(D.rewL_lr_blc{s,e,b}...
                                (bslvec)))./std(D.rewL_lr_blc{s,e,b}(bslvec));
                            tmp_df = tmp_h-tmp_l;
                        end

                        % Store for all delays and frequencies
                        delay(i).freq(b).etw_hr{dcnt,ecnt} = mean(tmp_h(ebinidx1:ebinidx2),'omitnan');
                        delay(i).freq(b).etw_lr{dcnt,ecnt} = mean(tmp_l(ebinidx1:ebinidx2),'omitnan');
                        delay(i).freq(b).etw_df{dcnt,ecnt} = mean(tmp_df(ebinidx1:ebinidx2),'omitnan');

                        delay(i).freq(b).ltw_hr{dcnt,ecnt} = mean(tmp_h(lbinidx1:lbinidx2),'omitnan');
                        delay(i).freq(b).ltw_lr{dcnt,ecnt} = mean(tmp_l(lbinidx1:lbinidx2),'omitnan');
                        delay(i).freq(b).ltw_df{dcnt,ecnt} = mean(tmp_df(lbinidx1:lbinidx2),'omitnan');
                    end
                    delay(i).freq(b).etw_hr{dcnt,1} = rname;
                    delay(i).freq(b).etw_lr{dcnt,1} = rname;
                    delay(i).freq(b).etw_df{dcnt,1} = rname;
                    delay(i).freq(b).ltw_hr{dcnt,1} = rname;
                    delay(i).freq(b).ltw_lr{dcnt,1} = rname;
                    delay(i).freq(b).ltw_df{dcnt,1} = rname;

                    delay(i).freq(b).etw_hr{dcnt,2} = s;
                    delay(i).freq(b).etw_lr{dcnt,2} = s;
                    delay(i).freq(b).etw_df{dcnt,2} = s;
                    delay(i).freq(b).ltw_hr{dcnt,2} = s;
                    delay(i).freq(b).ltw_lr{dcnt,2} = s;
                    delay(i).freq(b).ltw_df{dcnt,2} = s;

                    delay(i).freq(b).etw_hr{dcnt,3} = Delay_len;
                    delay(i).freq(b).etw_lr{dcnt,3} = Delay_len;
                    delay(i).freq(b).etw_df{dcnt,3} = Delay_len;
                    delay(i).freq(b).ltw_hr{dcnt,3} = Delay_len;
                    delay(i).freq(b).ltw_lr{dcnt,3} = Delay_len;
                    delay(i).freq(b).ltw_df{dcnt,3} = Delay_len;

                    numTr = D.behm.numHighTr + D.behm.numLowTr;
                    delay(i).freq(b).etw_hr{dcnt,4} = D.behm.numHighTr/numTr;
                    delay(i).freq(b).etw_lr{dcnt,4} = D.behm.numHighTr/numTr;
                    delay(i).freq(b).etw_df{dcnt,4} = D.behm.numHighTr/numTr;
                    delay(i).freq(b).ltw_hr{dcnt,4} = D.behm.numHighTr/numTr;
                    delay(i).freq(b).ltw_lr{dcnt,4} = D.behm.numHighTr/numTr;
                    delay(i).freq(b).ltw_df{dcnt,4} = D.behm.numHighTr/numTr;
                    dcnt=dcnt+1;
                end
            end
        end
    end
end
toc
disp('DD rewL plotting preprocessed');

%% Store Data as Table All Delay Version
sdir = 'E:\DD_PhysProcessed\ExcelStatsOutput\';
if ~exist(sdir,'dir')
    mkdir(sdir)
end
saveVars=1;
tic
tbltitle = [{'Animal'},{'Session'},{'DelayLength'},{'HR_choice'},{'M2'},{'A32D'},{'A32V'},{'vOFC'},{'LFC'}...
    {'ALM'},{'LOFC'},{'Ains'},{'VMS'},{'NAcC'},{'NAcS'},{'BLA'}];
% Convert cell arrays to tables 
for d = 1:length(delay)
    for b = 1:6
        etw_hr{d,b} =  cell2table(delay(d).freq(b).etw_hr);
        etw_lr{d,b} =  cell2table(delay(d).freq(b).etw_lr);
        etw_df{d,b} =  cell2table(delay(d).freq(b).etw_df);
        ltw_hr{d,b} =  cell2table(delay(d).freq(b).ltw_hr);
        ltw_lr{d,b} =  cell2table(delay(d).freq(b).ltw_lr);
        ltw_df{d,b} =  cell2table(delay(d).freq(b).ltw_df);

        % Change Table Title Names
        etw_hr{d,b}.Properties.VariableNames = tbltitle;
        etw_lr{d,b}.Properties.VariableNames = tbltitle;
        etw_df{d,b}.Properties.VariableNames = tbltitle;
        ltw_hr{d,b}.Properties.VariableNames = tbltitle;
        ltw_lr{d,b}.Properties.VariableNames = tbltitle;
        ltw_df{d,b}.Properties.VariableNames = tbltitle;
    end
end

dvec = [1,3,5,6,2,4]; % Order Delays are stored in structure
delay_list = ["0.5s","1s","2s","5s","10s","20s"];
if(saveVars)
    % Only store late time window(0-100ms) post reward window feedback
    for b = 1:6
        filename=char(strcat(freqs(b),'_allDelaysLTW_RewEs_stats.xlsx'));
        for d = 1:length(delay)
            dlen = ltw_hr{d,b}.('DelayLength');
            dlay = dlen(1);
            shname = strcat("Delay_",num2str(dlay),"s_LTW_HR");
            writetable(ltw_hr{d,b},fullfile(sdir,filename),'FileType','spreadsheet','Sheet',shname);
            shname = strcat("Delay_",num2str(dlay),"s__LTW_LR");
            writetable(ltw_lr{d,b},fullfile(sdir,filename),'FileType','spreadsheet','Sheet',shname);
            shname = strcat("Delay_",num2str(dlay),"s_LTW_DF");
            writetable(ltw_df{d,b},fullfile(sdir,filename),'FileType','spreadsheet','Sheet',shname);
            shname = strcat("Delay_",num2str(dlay),"s_ETW_HR");
            writetable(etw_hr{d,b},fullfile(sdir,filename),'FileType','spreadsheet','Sheet',shname);
            shname = strcat("Delay_",num2str(dlay),"s__ETW_LR");
            writetable(etw_lr{d,b},fullfile(sdir,filename),'FileType','spreadsheet','Sheet',shname);
            shname = strcat("Delay_",num2str(dlay),"s_ETW_DF");
            writetable(etw_df{d,b},fullfile(sdir,filename),'FileType','spreadsheet','Sheet',shname);
        end
        % Remove Extra Blank Sheets
        excelFileName=filename;
        excelFilePath = sdir;
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
toc
cd(sdir)
disp('DD Phys All delay LTW Stats Excel data stored');


%% Store Data as Table Grouped Delay Version
sdir = 'E:\DD_PhysProcessed\ExcelStatsOutput\';
if ~exist(sdir,'dir')
    mkdir(sdir)
end
saveVars=1;
tic
tbltitle = [{'Animal'},{'Session'},{'DelayLength'},{'HRchoice'},{'M2'},{'A32D'},{'A32V'},{'vOFC'},{'LFC'}...
    {'ALM'},{'LOFC'},{'Ains'},{'VMS'},{'NAcC'},{'NAcS'},{'BLA'}];
% Convert cell arrays to tables 
beta_ETW_ld_lr = cell2table(etwBeta_ld_lr);
beta_ETW_ld_hr = cell2table(etwBeta_ld_hr);
beta_ETW_ld_df = cell2table(etwBeta_ld_df);
beta_ETW_md_lr = cell2table(etwBeta_md_lr);
beta_ETW_md_hr = cell2table(etwBeta_md_hr);
beta_ETW_md_df = cell2table(etwBeta_md_df);
beta_ETW_hd_lr = cell2table(etwBeta_hd_lr);
beta_ETW_hd_hr = cell2table(etwBeta_hd_hr);
beta_ETW_hd_df = cell2table(etwBeta_hd_df);

beta_LTW_ld_lr = cell2table(ltwBeta_ld_lr);
beta_LTW_ld_hr = cell2table(ltwBeta_ld_hr);
beta_LTW_ld_df = cell2table(ltwBeta_ld_df);
beta_LTW_md_lr = cell2table(ltwBeta_md_lr);
beta_LTW_md_hr = cell2table(ltwBeta_md_hr);
beta_LTW_md_df = cell2table(ltwBeta_md_df);
beta_LTW_hd_lr = cell2table(ltwBeta_hd_lr);
beta_LTW_hd_hr = cell2table(ltwBeta_hd_hr);
beta_LTW_hd_df = cell2table(ltwBeta_hd_df);

% Change Table Title Names
beta_ETW_ld_lr.Properties.VariableNames = tbltitle;
beta_ETW_ld_hr.Properties.VariableNames = tbltitle;
beta_ETW_ld_df.Properties.VariableNames = tbltitle;
beta_ETW_md_lr.Properties.VariableNames = tbltitle;
beta_ETW_md_hr.Properties.VariableNames = tbltitle;
beta_ETW_md_df.Properties.VariableNames = tbltitle;
beta_ETW_hd_lr.Properties.VariableNames = tbltitle;
beta_ETW_hd_hr.Properties.VariableNames = tbltitle;
beta_ETW_hd_df.Properties.VariableNames = tbltitle;

beta_LTW_ld_lr.Properties.VariableNames = tbltitle;
beta_LTW_ld_hr.Properties.VariableNames = tbltitle;
beta_LTW_ld_df.Properties.VariableNames = tbltitle;
beta_LTW_md_lr.Properties.VariableNames = tbltitle;
beta_LTW_md_hr.Properties.VariableNames = tbltitle;
beta_LTW_md_df.Properties.VariableNames = tbltitle;
beta_LTW_hd_lr.Properties.VariableNames = tbltitle;
beta_LTW_hd_hr.Properties.VariableNames = tbltitle;
beta_LTW_hd_df.Properties.VariableNames = tbltitle;

if(saveVars)
    filename='Beta_LDMDHD_statsdataV2.xlsx';
    writetable(beta_ETW_ld_lr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_LD_LR');
    writetable(beta_ETW_ld_hr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_LD_HR');
    writetable(beta_ETW_ld_df,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_LD_DF');
    writetable(beta_ETW_md_lr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_MD_LR');
    writetable(beta_ETW_md_hr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_MD_HR');
    writetable(beta_ETW_md_df,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_MD_DF');
    writetable(beta_ETW_hd_lr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_HD_LR');
    writetable(beta_ETW_hd_hr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_HD_HR');
    writetable(beta_ETW_hd_df,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','ETW_HD_DF');
   
    writetable(beta_LTW_ld_lr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_LD_LR');
    writetable(beta_LTW_ld_hr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_LD_HR');
    writetable(beta_LTW_ld_df,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_LD_DF');
    writetable(beta_LTW_md_lr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_MD_LR');
    writetable(beta_LTW_md_hr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_MD_HR');
    writetable(beta_LTW_md_df,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_MD_DF');
    writetable(beta_LTW_hd_lr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_HD_LR');
    writetable(beta_LTW_hd_hr,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_HD_HR');
    writetable(beta_LTW_hd_df,fullfile(sdir,filename),'FileType','spreadsheet','Sheet','LTW_HD_DF');
    
    % Remove Extra Blank Sheets
    excelFileName=filename;
    excelFilePath = sdir;
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
toc
cd(sdir)
disp('DD Phys LD,MD,HD lrvshr Stats Excel data stored');

