% Delay Discounting Behavioral Data Excel Extraction
% For computatational modeling
% Extract behavior from all sessions with physiology
clearvars; close all; clc;
%% Store data across all trials in each session for each delay and each rat

rdir='E:\DD_PhysProcessed\rewL_processedData_Umesh';
sdir ='E:\DD_PhysProcessed\ExcelOutput_trialbytrial\'; % save into new directory
if~exist(sdir,'dir')
    mkdir(sdir);
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
%%
% switch based on which is running
tic
for i = 1:length(ddir) % Iterate through all delay lengths
    disp(ddir(i).name)
    cd([rdir '\' ddir(i).name]);
    adir = dir;
    adir(1:2) = [];
    delayL(i) = str2num(extractBefore(ddir(i).name,'ms')); %#ok<ST2NM>
    sddir = strcat(sdir,ddir(i).name,'\');
    if(~exist(sddir,'dir'))
        mkdir(sddir);
    end

    for j = 1:length(adir) %Iterate through all sessions(across all anis)
        an = adir(j).name;
        if(startsWith(an,'R'))
            rat = an;
            if ancnt>1
                clear D;
                rname{ancnt} = rat;
                ancnt=ancnt+1;
                %end
            else
                rname{ancnt} = rat;
                ancnt=ancnt+1;
            end
            %for kk=1:size(D.rewL_tr)

            cd([adir(j).folder '\' an]);
            andata=dir;
            andata(1:2)=[];
            for k=1:length(andata)
                if (contains(andata(k).name,'new'))
                    load(andata(k).name);
                    sescnt=0;
                    for l = 1:length(D.rewL_tr)
                        sescnt = l;
                        if (size(D.rewL_tr{l}{1},1))==32
                        D.An_phys{sescnt,1}=rat;
                        D.An_phys{sescnt,2}=delayL(i);
                        % date = extractBetween(an,'_DD_','_');
                        numTr = length(D.behm{sescnt}.resp);
                        trial = 1:numTr;
                        crat = rname{ancnt-1};
                        Subject = cell(numTr,1);
                        Subject(:) = {crat};
                        DelayL = cell(numTr,1);
                        DelayL(:)={D.behm{sescnt}.trType};

                        Session = cell(numTr,1);
                        Session(:) = {sescnt};
                        HRLatency = zeros(numTr,1);
                        LRLatency = zeros(numTr,1);
                        LRchoice = zeros(numTr,1);
                        HRchoice = zeros(numTr,1);
                        res_lr = D.behm{sescnt}.res_2;
                        res_hr = D.behm{sescnt}.res_4;
                        lrempt = cellfun(@isempty, res_lr);
                        hrempt = cellfun(@isempty, res_hr);
                        res_lr(lrempt) = {0};
                        res_hr(hrempt) = {0};
                        %                     LRLatency = cell2mat(res_lr);
                        %                     HRLatency = cell2mat(res_hr);
                        for tr = trial
                            np2v(tr) = ~isempty(D.behm{sescnt}.res_2{tr});
                            np4v(tr) = ~isempty(D.behm{sescnt}.res_4{tr});
                            if(np2v(tr)==1) % if NP2 is not empty
                                %                         tmp = D.behm.res_2{tr}-D.behm.;
                                %                            LRLatency(tr) = tmp(1);
                                LRchoice(tr) = 1;
                            elseif(np4v(tr)==1) % if NP4 is not empty
                                %                          tmp = D.behm.res_4{tr};
                                %                           HRLatency(tr) = tmp(1);
                                HRchoice(tr) = 1;
                            end
                            ChoiceLatency(tr,1)=D.behm{sescnt}.startRT(tr);
                            RewardLatency(tr,1)=D.behm{sescnt}.rewRT(tr);
                        end
                        filename = strcat(crat,'_delay_',num2str(delayL),'_',num2str(sescnt),'.xlsx');
                        TrialNum = trial';

                        dataT = table(Subject,DelayL,Session,TrialNum,...
                            LRchoice,HRchoice,ChoiceLatency,RewardLatency);
                        for aa = 1:32
                            Varnames{aa,1}=strcat('Elec',num2str(aa));

                        end

                        %delta
                        writetable(dataT,fullfile(sddir,char(filename)),'Sheet','Delta');
                        Elec=squeeze(nanmean(D.rewL_tr{sescnt}{1}(:,84:117,:),2))';
                        dataLFP = array2table(Elec);
                        dataLFP.Properties.VariableNames = Varnames;
                        writetable(dataLFP,fullfile(sddir,char(filename)),'Sheet','Delta','Range','I1');

                        %theta
                        writetable(dataT,fullfile(sddir,char(filename)),'Sheet','Theta');
                        Elec=squeeze(nanmean(D.rewL_tr{sescnt}{2}(:,84:117,:),2))';
                        dataLFP = array2table(Elec);
                        dataLFP.Properties.VariableNames = Varnames;
                        writetable(dataLFP,fullfile(sddir,char(filename)),'Sheet','Theta','Range','I1');

                        %alpha
                        writetable(dataT,fullfile(sddir,char(filename)),'Sheet','Alpha');
                        Elec=squeeze(nanmean(D.rewL_tr{sescnt}{3}(:,84:117,:),2))';
                        dataLFP = array2table(Elec);
                        dataLFP.Properties.VariableNames = Varnames;
                        writetable(dataLFP,fullfile(sddir,char(filename)),'Sheet','Alpha','Range','I1');

                        %beta
                        writetable(dataT,fullfile(sddir,char(filename)),'Sheet','Beta');
                        Elec=squeeze(nanmean(D.rewL_tr{sescnt}{4}(:,84:117,:),2))';
                        dataLFP = array2table(Elec);
                        dataLFP.Properties.VariableNames = Varnames;
                        writetable(dataLFP,fullfile(sddir,char(filename)),'Sheet','Beta','Range','I1');

                        %gamma
                        writetable(dataT,fullfile(sddir,char(filename)),'Sheet','Gamma');
                        Elec=squeeze(nanmean(D.rewL_tr{sescnt}{5}(:,84:117,:),2))';
                        dataLFP = array2table(Elec);
                        dataLFP.Properties.VariableNames = Varnames;
                        writetable(dataLFP,fullfile(sddir,char(filename)),'Sheet','Gamma','Range','I1');

                        %highgamma
                        writetable(dataT,fullfile(sddir,char(filename)),'Sheet','HighGamma');
                        Elec=squeeze(nanmean(D.rewL_tr{sescnt}{6}(:,84:117,:),2))';
                        dataLFP = array2table(Elec);
                        dataLFP.Properties.VariableNames = Varnames;
                        writetable(dataLFP,fullfile(sddir,char(filename)),'Sheet','HighGamma','Range','I1');

                        clear Subject DelayL Elec dataLFP Session TrialNum HRchoice LRchoice
                        clear ChoiceLatency RewardLatency behm
                        end
                    end
                    clear  D
                end
            end
        end
    end
end



cd(sdir)
disp('Done DD excel output across all trials(phys dates only)');

