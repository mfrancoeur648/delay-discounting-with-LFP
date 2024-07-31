% Delay Discounting Connectivity Script
% Calculate granger and WPLI connectivity across Delay conditions

clearvars;
close all;
%% Load 3Dim epoched Data
rdir='E:\dd_processed\trType_labeled\';
cd(rdir)
ddir = dir;
ddir(1:2)=[];

sdir = 'E:\DD_PhysProcessed\Connectivity\';
if~exist(sdir,'dir')
   mkdir(sdir);
end 
pdir1 = 'E:\DD_PhysProcessed\Connectivity\GrangerPlots\';
if ~exist(pdir1,'dir')
    mkdir(pdir1)
end

load('E:\DD_PhysProcessed\times_freqs.mat');
% Set Desired Time Window for plotting 

tw1 = -500;
tw2 = 2500;

twa1 = tw1-30;
twa2 = tw2+30;
tt = t(t>twa1 & t<twa2);
ttvec = find(t==min(tt)):find(t==max(tt));

% Early Time Window (pre response)
etw1 = -2000;
etw2 = -1000;
etwa1 = etw1-30;
etwa2 = etw2+30;
tt_e = t(t>etwa1 & t<etwa2);
etvec = find(t==min(tt_e)):find(t==max(tt_e));
ebinidx1 = min(etvec);
ebinidx2 = max(etvec);

% Middle Time Window
mtw1 = 0;
mtw2 = 1000;
mtwa1 = mtw1-30;
mtwa2 = mtw2+30;
tt_m = t(t>mtwa1 & t<mtwa2);
mtvec = find(t==min(tt_m)):find(t==max(tt_m));
mbinidx1 = min(mtvec);
mbinidx2 = max(mtvec);

% Late Time Window
ltw1 = 1000;
ltw2 = 2000;
ltwa1 = ltw1-30;
ltwa2 = ltw2+30;
tt_l = t(t>ltwa1 & t<ltwa2);
ltvec = find(t==min(tt_l)):find(t==max(tt_l));
lbinidx1 = min(ltvec);
lbinidx2 = max(ltvec);

% Scale to -3000-4000 time vector (behLFP)

etw1 = etw1+3000;
etw2 = etw2+3000;

mtw1 = mtw1+3000;
mtw2 = mtw2+3000;

ltw1 = ltw1+3000;
ltw2 = ltw2+3000;

saveVars=0;
saveFigs=1;
plotFigs=1;
eVec = [28,27,26,25,31,32,29,30,6,5,1,20]; 
runConnectivity=0;
% switch based on which is running    
ancnt=1; %reset for each delay length 
sescnt=0;
tic
for i = 1:length(ddir) % Iterate through all delay lengths 
    disp(ddir(i).name)         
    cd([rdir '\' ddir(i).name]);
    adir = dir;
    adir(1:2) = []; 
    delayL(i) = str2num(extractBefore(ddir(i).name,'msH')); %#ok<ST2NM>
    Delay(i).length = delayL(i);

    for j = 1:length(adir) %Iterate through all sessions(across all anis)
        an = adir(j).name;
        if(startsWith(an,'R')&&contains(an,'_DD_'))                
            rat = extractBefore(an,'_DD_');
            if ancnt>1
                if ~strcmp(rat,rname{ancnt-1})
                    if(saveVars)
                        if(exist('D','var'))
%                             savedir = strcat(sdir,num2str(delayL(i)),'msDelay\',rname{ancnt-1},'\');
%                             if~exist(savedir,'dir')
%                                mkdir(savedir);
%                             end 
%                             fname = strcat(rname{ancnt-1},'_rewLbinned.mat');
%                             save(fullfile(char(savedir),char(fname)),'D','-v7.3');
                            disp(rname{ancnt-1});
                        end
                    end
                    clear D;
                    rname{ancnt} = rat;
                    ancnt=ancnt+1;
                    sescnt=0;
                end
            else
                rname{ancnt} = rat;
                ancnt=ancnt+1;
            end
                
            cd([adir(j).folder '\' an]);
            if(exist('epochdata.mat','file') && ....
                    exist('extracted_data.mat','file'))   
                sescnt = sescnt+1;
                D.An_phys{sescnt,1}=rat;
                D.An_phys{sescnt,2}=delayL(i);
                load('extracted_data.mat');
                date = extractBetween(an,'_DD_','_');
                % store behavioral data across sessions
                numTr = behm.numHighTr + behm.numLowTr;
                Del(i).hr_ch{ancnt,sescnt} = behm.numHighTr/numTr;
                if(runConnectivity)                
                    load('epochdata.mat');
                    E = size(behLFP.rew,1);
                    L = size(behLFP.rew,3);
                    lb = [0,4,8,13,30,70];
                    hb = [4,8,13,30,50,120];
                    for e = 1:E % iterate through all electrodes
                        lrcnt = 0;
                        hrcnt = 0;
                        % trial type sort raw epoched data across TWs
                        for tr = 1:L % 
                           % Store low and high reward vectors
                            % across all trials 
                            np2v(tr) = ~isempty(behm.res_2{tr});
                            np4v(tr) = ~isempty(behm.res_4{tr});
                            if(np2v(tr)==1) % if NP2 is not empty 
                                lrcnt = lrcnt+1;
                                TW(1).lr(e,:,lrcnt) = behLFP.rew(e,etw1:etw2,tr); 
                                TW(2).lr(e,:,lrcnt) = behLFP.rew(e,mtw1:mtw2,tr); 
                                TW(3).lr(e,:,lrcnt) = behLFP.rew(e,ltw1:ltw2,tr); 
                            elseif(np4v(tr)==1) % if NP4 is not empty
                                hrcnt = hrcnt+1;
                                TW(1).hr(e,:,hrcnt) = behLFP.rew(e,etw1:etw2,tr); 
                                TW(2).hr(e,:,hrcnt) = behLFP.rew(e,mtw1:mtw2,tr); 
                                TW(3).hr(e,:,hrcnt) = behLFP.rew(e,ltw1:ltw2,tr);
                            end
                        end
                    end
                    % Store separate struct for 12 rew related elecs
                    ecnt=0;
                    for e = eVec 
                        ecnt= ecnt+1;
                        lrcnt = 0;
                        hrcnt = 0;
                        % trial type sort raw epoched data across TWs
                        for tr = 1:L % 
                           % Store low and high reward vectors
                            % across all trials 
                            np2v(tr) = ~isempty(behm.res_2{tr});
                            np4v(tr) = ~isempty(behm.res_4{tr});
                            if(np2v(tr)==1) % if NP2 is not empty 
                                lrcnt = lrcnt+1;
                                rewTW(1).lr(ecnt,:,lrcnt) = behLFP.rew(e,etw1:etw2,tr); 
                                rewTW(2).lr(ecnt,:,lrcnt) = behLFP.rew(e,mtw1:mtw2,tr); 
                                rewTW(3).lr(ecnt,:,lrcnt) = behLFP.rew(e,ltw1:ltw2,tr); 
                            elseif(np4v(tr)==1) % if NP4 is not empty
                                hrcnt = hrcnt+1;
                                rewTW(1).hr(ecnt,:,hrcnt) = behLFP.rew(e,etw1:etw2,tr); 
                                rewTW(2).hr(ecnt,:,hrcnt) = behLFP.rew(e,mtw1:mtw2,tr); 
                                rewTW(3).hr(ecnt,:,hrcnt) = behLFP.rew(e,ltw1:ltw2,tr);
                            end
                        end
                    end

                    % Connectivity Calculation for each session
                    for x = 1:length(TW) % iterate though all TWs
                        trdat{1} = TW(x).lr;
                        trdat{2} = TW(x).hr;

                        rewdat{1} = rewTW(x).lr;
                        rewdat{2} = rewTW(x).hr;

                        for d = 1:length(trdat)  % iterate through all trial types

                            % remove nan data - missing electrodes
                            if any(isnan(trdat{d}(:,1,1)))
                                badE = find(isnan(trdat{d}(:,1,1)));
                                for bd = badE
                                    trdat{d}(bd,:,:) = [];
                                end
                            end
                            if any(isnan(rewdat{d}(:,1,1)))
                                badE2 = find(isnan(rewdat{d}(:,1,1)));
                                for bd = badE2
                                    rewdat{d}(bd,:,:) = [];
                                end
                            end
                            % Computation of multivariate autoregressive
                            % model
                            if size(trdat{d},1)<20
                                continue;
                            end
                            numtp = size(trdat{d},2);
                            numtp1 = size(rewdat{d},2);

                            sr=1000;
                            times=(0.001:numtp)/sr;
                            times1 = (0.001:numtp1)/sr;
    %                             times = 1:size(trdat{d},2);
                            data = [];
                            data.ntrials = size(trdat{d},3);
                            data.nsignal = size(trdat{d},1); % electrodes
                            data.fsample = 1000; % check
                            data.triallength = 1;

                            data1 = [];
                            data1.ntrials = size(rewdat{d},3);
                            data1.nsignal = size(rewdat{d},1); % electrodes
                            data1.fsample = 1000; % check
                            data1.triallength = 1;

                            for k = 1:size(trdat{d},3)
                                data.trial{1,k} = double(trdat{d}(:,:,k));
                                data.time{1,k} = double(times);
                            end
                            for k = 1:size(rewdat{d},3)
                                data1.trial{1,k} = double(rewdat{d}(:,:,k));
                                data1.time{1,k} = double(times1);
                            end
                            vcnt=1;
                            for v = 1:32
                                if exist('badE','var')
                                    if ismember(v,badE)
                                        continue;
                                    end
                                end
                                data.label{vcnt,1}=strcat('Electrode',num2str(v));
                                vcnt = vcnt+1;
                            end
                            vcnt = 1;
                            for v = 1:12
                                if exist('badE2','var')
                                    if ismember(v,badE2)
                                        continue;
                                    end
                                end
                                data1.label{vcnt,1} = strcat('Electrode',num2str(eVec(v)));
                                vcnt = vcnt+1;
                            end
    %                             cfg         = [];
    %                             cfg.order   = 10;
    %                             cfg.toolbox = 'bsmart';
    %                             cfg.keeptrials = 'yes';
    %                             %keyboard;
    %                             mdata       = ft_mvaranalysis(cfg, data);
    %                             mdata1      = ft_mvaranalysis(cfg, data1);
    %                             mdata.label = data.label;
    %                             mdata1.label = data1.label;
    %                             %keyboard;
                            % do spectral analysis 
    %                             cfg = [];
    %                             cfg.method    = 'mtmfft';
    %                             cfg.output    = 'fourier';
    %                             cfg.foilim    = [0 200];
    %                             cfg.tapsmofrq = 5;
    %                             freq          = ft_freqanalysis(cfg, data{d});
    %                             fd            = ft_freqdescriptives(cfg, freq);
    %                             figure;
    %                             plot(fd.freq, fd.powspctrm);
    %                             set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
    %                             title('power spectrum');

                            % Non parametric computation of the
                            % cross-spectral density matrix
                            cfg           = [];
                            cfg.method    = 'mtmfft';
                            cfg.taper     = 'dpss';
                            cfg.output    = 'fourier';
                            cfg.tapsmofrq = 2;
                            try
                                freq          = ft_freqanalysis(cfg, data);
                                freq1         = ft_freqanalysis(cfg, data1);
                            catch ME
                                if(strcmp(ME,"data = ft_checkdata(data, ...'datatype', {'raw', 'raw+comp', 'mvar'},'feedback', cfg.feedback, 'hassampleinfo', 'yes');"))
                                    continue;
                                end
                            end
                            fd            = ft_freqdescriptives(cfg, freq);
                            figure;
                            plot(fd.freq, fd.powspctrm);
                            set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
                            title('power spectrum');
    %                             cfg        = [];
    %                             cfg.method = 'mvar';
    %                             mfreq      = ft_freqanalysis(cfg, mdata);

                            % Computation and visualization of coherence
                            % coefficient
                            cfg           = [];
                            cfg.method    = 'coh';
                            coh           = ft_connectivityanalysis(cfg, freq);
    %                             cohm          = ft_connectivityanalysis(cfg, mfreq);

                            cfg           = [];
                            cfg.parameter = 'cohspctrm';
                            cfg.zlim      = [0 1];
                            ft_connectivityplot(cfg, coh);

                            % Computation and visualization of granger
                            % connectiity
                            cfg           = [];
                            cfg.method    = 'granger';
                            % save granger for both 32 region and 12
                            % regions across all delays
                            Delay(i).granger{ancnt,sescnt,x,d}= ft_connectivityanalysis(cfg, freq);
                            Delay(i).granger_rew{ancnt,sescnt,x,d}= ft_connectivityanalysis(cfg, freq1);

                            cfg           = [];
                            cfg.parameter = 'grangerspctrm';
                            cfg.zlim      = [0 1];
                            ft_connectivityplot(cfg, Delay(i).granger{ancnt,sescnt,x,d});
                            flabel = strcat('Rat',{' '},rname{ancnt-1},'Session',{' '},...
                                num2str(sescnt),'All Elec Granger');
                            if(saveFigs)
                                saveas(gcf,fullfile(pdir1,char(flabel)),'fig');
                            end
                            ft_connectivityplot(cfg, Delay(i).granger_rew{ancnt,sescnt,x,d});
                            flabel = strcat('Rat',{' '},rname{ancnt-1},'Session',{' '},...
                                num2str(sescnt),'Rew Related Granger');
                            if(saveFigs)
                                saveas(gcf,fullfile(pdir1,char(flabel)),'fig');
                            end

                            % Computation and Visualization of WPLI
                            cfg = [];
                            cfg.method = 'wpli_debiased';
                            % save wpli fpor both 32 region and 12 rew regions
                            % across all delays
                            Delay(i).wpli{ancnt,sescnt,x,d} = ft_connectivityanalysis(cfg, freq);
                            Delay(i).wpli_rew{ancnt,sescnt,x,d} = ft_connectivityanalysis(cfg, freq1);

    %                           wpli = ft_checkdata(wpli, 'cmbrepresentation', 'full','datatype','freq');
    %                             ft_connectivityplot(cfg, wpli);
                            clear badE badE2
                            close all;
                        end
                    end
                end
            end
        end
    end
end
if(saveVars)
    % if(exist('Delay','var'))
      
     %   fname = 'DD_granger_wpli2.mat';
     %   save(fullfile(sdir,fname),'Delay','Del','-v7.3');
    %end
    
    if(exist('Del','var'))
        fname = 'DD_granger_wpli2.mat';
        save(fullfile(sdir,fname),'Del','-v7.3');
    end
end
toc
disp('Done processing granger and WPLI');      


%%
clearvars -except Del
ldir = 'E:\DD_PhysProcessed\Connectivity';
tic
load(fullfile(ldir,'DD_granger_wpli.mat'));
load(fullfile(ldir,'DD_granger_wpli2.mat'));
toc

%% Averaging Connectivity output across all sessions and animals
% Store data concatenated across all session/animals pre avg
% Extract out connectivity across frequency bands
for i = 1:length(Delay)
    % Average across all sesssions scross all animals 

    freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];

    elecs = ["NAcSh","A33","A24a","A24b","NAcCr","VMS","DMS","M1p","V1","V1","DorSub",...
        "DG","PPCx","CA1","CA3","STN","DLS","DLS","CE Am","BLA","A30c","A29c","MD Th",...
        "CE Th","vOFC","A32V","A32D","M2","L OFC","vAI","LFC","M1a"];
    A = size(Delay(i).granger,1); % animals
    S = size(Delay(i).granger,2); % sessions
    numTWs = size(Delay(i).granger,3); %time windows (etw, mtw, ltw)
    numTrTypes = size(Delay(i).granger,4); % tr types we are looking at
    lb = [0,4,8,13,30,70];
    hb = [4,8,13,30,50,120];
    for tw = 1:numTWs
        for tt = 1:numTrTypes
            r_scnt =0; % avergage across sessions and animals
            ae_scnt=0;
            for a = 1:A
                for s = 1:S 
                    % take mean within granger matrix output
                    % Avg across granger spectrum and store elec matrix
                    
                    if ~isempty(Delay(i).granger{a,s,tw,tt})
                        if(size(Delay(i).granger{a,s,tw,tt}.grangerspctrm,1)>31) 
                            % dont include sessions with less than 32 elecs
                            ae_scnt = ae_scnt+1;
                            for b = 1:numel(lb)
                                bvec = find(Delay(i).granger{a,s,tw,tt}.freq>lb(b) & ...
                                    Delay(i).granger{a,s,tw,tt}.freq<hb(b));
                                Delay(i).F(b).TW(tw).TrType(tt).gtmp(ae_scnt,:,:) = ...
                                    mean(Delay(i).granger{a,s,tw,tt}.grangerspctrm(:,:,bvec),3); 
                                Delay(i).F(b).TW(tw).TrType(tt).wpltmp(ae_scnt,:,:) = ...
                                    mean(Delay(i).wpli{a,s,tw,tt}.wpli_debiasedspctrm(:,:,bvec),3); 
                            end
                            if(~isempty(Del(i).hr_ch{a,s}))
                                Delay(i).hr_ch_all{ae_scnt} = Del(i).hr_ch{a,s};
                            end
                        end
                    end
                    if ~isempty(Delay(i).granger_rew{a,s,tw,tt})
                        if(size(Delay(i).granger_rew{a,s,tw,tt}.grangerspctrm,1)>11) 
                            % dont include sessions with less than 12 elecs
                            r_scnt = r_scnt+1;
                            for b = 1:numel(lb)
                                bvec = find(Delay(i).granger_rew{a,s,tw,tt}.freq>lb(b) & ...
                                        Delay(i).granger_rew{a,s,tw,tt}.freq<hb(b));
                                Delay(i).F(b).TW(tw).TrType(tt).grtmp(r_scnt,:,:) = ...
                                    mean(Delay(i).granger_rew{a,s,tw,tt}.grangerspctrm(:,:,bvec),3);
                                Delay(i).F(b).TW(tw).TrType(tt).wplrtmp(r_scnt,:,:) = ...
                                    mean(Delay(i).wpli_rew{a,s,tw,tt}.wpli_debiasedspctrm(:,:,bvec),3);
                            end
                            if(~isempty(Del(i).hr_ch{a,s}))
                                Delay(i).hr_ch_allR{r_scnt} = Del(i).hr_ch{a,s};
                            end
                        end
                    end
                end
            end
            for b =1:numel(lb)
                Delay(i).F(b).gmat{tw,tt} = ...
                    squeeze(mean(Delay(i).F(b).TW(tw).TrType(tt).gtmp,1,'omitnan'));
                Delay(i).F(b).grmat{tw,tt} = ...
                    squeeze(mean(Delay(i).F(b).TW(tw).TrType(tt).grtmp,1,'omitnan'));
                Delay(i).F(b).wplmat{tw,tt} = ...
                    squeeze(mean(Delay(i).F(b).TW(tw).TrType(tt).wpltmp,1,'omitnan'));
                Delay(i).F(b).wplrmat{tw,tt} = ...
                    squeeze(mean(Delay(i).F(b).TW(tw).TrType(tt).wplrtmp,1,'omitnan'));
            end
        end
    end
end
disp('DD Granger and WPLI averaged across all animals and sessions');

%% WPLI Connecitivty plots 

pdir = 'E:\DD_PhysProcessed\Connectivity\Results\GrangerPlots\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
pdir2 = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\';
if ~exist(pdir2,'dir')
    mkdir(pdir2)
end
saveFigs =1;
elecs = ["NAcSh","A33","A24a","A24b","NAcCr","VMS","DMS","M1p","V1","V1","DorSub",...
    "DG","PPCx","CA1","CA3","STN","DLS","DLS","CE Am","BLA","A30c","A29c","MD Th",...
    "CE Th","vOFC","A32V","A32D","M2","L OFC","vAI","LFC","M1a"];
rew_elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];

trTypes = ["LR","HR"];
twlab = ["ETW","MTW","LTW"];
for i = 1%:length(Delay)
    for b = 4
        for tt = 1:2% iterate hr and lr 
            % 32 elec granger plotting
            %tmp = (gmat{tw,tt}-min(gmat{tw,tt}))/(max(gmat{tw,tt})-min(gmat{tw,tt}));
    %             schemaball(elecs,Delay(i).F(b).gmat{tw,tt})
    %             flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Schemaball');
    %             title(flabel);
    %             if(saveFigs)
    %                 saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    % %             end
    %         figure;
    %         subplot(221);
    %         imagesc(Delay(i).F(b).gmat{3,tt});colorbar; caxis([-0.5,0.5]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('Late TW');
    %         subplot(222);
    %         imagesc(Delay(i).F(b).gmat{2,tt});colorbar; caxis([-0.5,0.5]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('Mid TW');
    %         subplot(223)
    %         imagesc(Delay(i).F(b).gmat{3,tt}-Delay(i).F(b).gmat{1,tt});colorbar; caxis([-0.5,0.5]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr LTW');
    %         subplot(224)
    %         imagesc(Delay(i).F(b).gmat{2,tt}-Delay(i).F(b).gmat{1,tt}); colorbar; caxis([-0.5,0.5]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr MTW');
    %         
    %         figure;
    %         subplot(221);
    %         imagesc(Delay(i).F(b).grmat{3,tt});colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('Late TW');
    %         subplot(222);
    %         imagesc(Delay(i).F(b).grmat{2,tt});colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('Mid TW');
    %         subplot(223)
    %         imagesc(Delay(i).F(b).grmat{3,tt}-Delay(i).F(b).grmat{1,tt});colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr LTW');
    %         subplot(224)
    %         imagesc(Delay(i).F(b).grmat{2,tt}-Delay(i).F(b).grmat{1,tt}); colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr MTW');
    %         
            figure;
            subplot(221);
            imagesc(Delay(i).F(b).wplrmat{3,tt});colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Late TW(1 to 2s)');
            subplot(222);
            imagesc(Delay(i).F(b).wplrmat{2,tt});colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Mid TW(0 to 1s)');
            subplot(223)
            imagesc(Delay(i).F(b).wplrmat{3,tt}-Delay(i).F(b).wplrmat{1,tt});colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Baseline Corrected LTW');
            subplot(224)
            imagesc(Delay(i).F(b).wplrmat{2,tt}-Delay(i).F(b).wplrmat{1,tt}); colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Baseline Corrected MTW');
            flabel = strcat(num2str(Delay(i).length),'s',{' '},'Delay',{' '},freqs(b),...
                {' '},trTypes(tt),{' '},'Rew Regions WPLI Connectivity');
            %suptitle(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
            end
%            close;
    %         
    %         figure;
    %         subplot(221);
    %         imagesc(Delay(i).F(b).wplmat{3,tt});colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('Late TW');
    %         subplot(222);
    %         imagesc(Delay(i).F(b).wplmat{2,tt});colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel', elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel', elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('Mid TW');
    %         subplot(223)
    %         imagesc(Delay(i).F(b).wplmat{3,tt}-Delay(i).F(b).wplmat{1,tt});colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr LTW');
    %         subplot(224)
    %         imagesc(Delay(i).F(b).wplmat{2,tt}-Delay(i).F(b).wplmat{1,tt}); colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr MTW');
    %         suptitle('WPLI');

    %         subplot(235)
    %         imagesc(Delay(i).F(b).wplmat{1,tt}); colorbar; caxis([-0.25,0.25]);
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    % %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Connectivty');
    %         title('BSL Corr MTW');
    %         suptitle('WPLI');
    %         
    %         



    %         if(saveFigs)
    %             saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    %         end
    %         % 12 rew elec granger plotting
    % %             figure;
    % %             schemaball(rew_elecs,Delay(i).F(b).grmat{tw,tt})
    % %             flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Rew Elecs Schemaball');
    % %             title(flabel);
    % %             if(saveFigs)
    % %                 saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    % %             end
    %         figure;
    %         imagesc(grmat{tw,tt})
    %         set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
    %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'Granger Rew Elecs Connectivty');
    %         title(flabel);
    %         if(saveFigs)
    %             saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    %         end
    % 
    %         % WPLI plotting
    %         % 32 elec wpli plotting
    %         figure;
    %         schemaball(elecs,Delay(i).F(b).wplmat{tw,tt})
    %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'WPLI Schemaball');
    %         title(flabel);
    %         if(saveFigs)
    %             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
    %         end
    %         figure;
    %         imagesc(Delay(i).F(b).wplmat{tw,tt})
    %         set(gca, 'XTick',1:numel(elecs), 'XTickLabel',elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(elecs), 'YTickLabel',elecs, 'YTickLabelRotation',45);
    %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'WPLI Connectivty');
    %         title(flabel);
    %         if(saveFigs)
    %             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
    %         end
    %         % 12 rew elec granger plotting
    %         figure;
    %         schemaball(rew_elecs,Delay(i).F(b).wplrmat{tw,tt})
    %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'WPLI Rew Elecs Schemaball');
    %         title(flabel);
    %         if(saveFigs)
    %             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
    %         end
    %         figure;
    %         imagesc(Delay(i).F(b).wplrmat{tw,tt})
    %         set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
    %         set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
    %         flabel = strcat(twlab(tw),{' '},trTypes(tt),{' '},'WPLI Rew Elecs Connectivty');
    %         title(flabel);
    %         if(saveFigs)
    %             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
    %         end
    %         pause;
        end
    end
end
disp('Done plotting connectivity');
cd(pdir)
%% Connectivity Stats mean across sessions
for i = 1:length(Delay)
    for b =1:6
        for tt = 1:2
            for tw =1:3
                [~,Delay(i).F(b).TW(tw).TrType(tt).wpl_tt] = ttest(Delay(i).F(b).TW(tw).TrType(tt).wpltmp);
                [~,Delay(i).F(b).TW(tw).TrType(tt).wplr_tt] = ttest(Delay(i).F(b).TW(tw).TrType(tt).wplrtmp);
            end
        
        bsldf=Delay(i).F(b).TW(3).TrType(tt).wplrtmp-Delay(i).F(b).TW(1).TrType(tt).wplrtmp;
        Delay(i).F(b).TrType(tt).wplr_ltw_bsl = mean(bsldf,1,'omitnan');
        [~,Delay(i).F(b).TrType(tt).wplr_ltw_bsl_p] = ttest(bsldf);
        
        bsldf2=Delay(i).F(b).TW(2).TrType(tt).wplrtmp-Delay(i).F(b).TW(1).TrType(tt).wplrtmp;
        Delay(i).F(b).TrType(tt).wplr_mtw_bsl = mean(bsldf2,1,'omitnan');
        [~,Delay(i).F(b).TrType(tt).wplr_mtw_bsl_p] = ttest(bsldf2);
        
        % Paired ttests HR-LR 
%         Delay(i).F(b).TW(tw).wpl_hrlr_tt = squeeze(ttest(Delay(i).F(b).TW(tw).TrType(1).wplrtmp,...
%            Delay(i).F(b).TW(tw).TrType(2).wplrtmp));
%         Delay(i).F(b).TW(tw).wpl_hrlr_tt = squeeze(ttest(Delay(i).F(b).TW(tw).TrType(1).wplrtmp,...
%             Delay(i).F(b).TW(tw).TrType(2).wplrtmp));
%         
%         % Wilcoxon Rank Sum test
%         Delay(i).F(b).TW(tw).wpl_hrlr_wrs = squeeze(ranksum(Delay(i).F(b).TW(tw).TrType(1).wplrtmp,...
%            Delay(i).F(b).TW(tw).TrType(2).wplrtmp));
%         Delay(i).F(b).TW(tw).wpl_hrlr_wrs = squeeze(ranksum(Delay(i).F(b).TW(tw).TrType(1).wplrtmp,...
%             Delay(i).F(b).TW(tw).TrType(2).wplrtmp));
%                 
        end
        LL = length(Delay(i).F(b).TW(2).TrType(1).wplrtmp);
        % HR-LR trial type differerences (MTW - 0-1s DD important window)
        deldif_r = Delay(i).F(b).TW(2).TrType(2).wplrtmp(1:LL,:,:)-Delay(i).F(b).TW(2).TrType(1).wplrtmp;
        Delay(i).F(b).wplr_hrlrdf = mean(deldif_r,1,'omitnan');
        [~,Delay(i).F(b).wplr_hrlrdf_p] = ttest(deldif_r);
        
%         % HR-LR trial type differerences (MTW - 0-1s DD important window)
%         deldif = Delay(i).F(b).TW(2).TrType(2).wpltmp-Delay(i).F(b).TW(2).TrType(1).wpltmp;
%         Delay(i).F(b).wpl_hrlrdf = mean(deldif,1,'omitnan');
%         [~,Delay(i).F(b).wpl_hrlrdf_p] = ttest(deldif);
    end
end
disp('Performed DD Connectivity Stats across Sessions');

keyboard;
%% Thresholded Connectivity plotting 
                            
pdir2 = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\Thresholded\';
if ~exist(pdir2,'dir')
    mkdir(pdir2) 
end
for i = 1:length(Delay)
    for b = 4
        for tt = 1:2 % iterate hr and lr
            figure;
            subplot(221);
            tmp1 = Delay(i).F(b).wplrmat{3,tt};
            tmp1(fdr(Delay(i).F(b).TW(3).TrType(tt).wplr_tt)>0.05) = 0;
            imagesc(tmp1);colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Late TW(1 to 2s)');
            tmp2 = Delay(i).F(b).wplrmat{2,tt};
            tmp2(fdr(Delay(i).F(b).TW(2).TrType(tt).wplr_tt)>0.05) = 0;
            subplot(222);
            imagesc(tmp2);colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Mid TW(0 to 1s)');
            
            subplot(223)
            tmp3 = Delay(i).F(b).TrType(tt).wplr_ltw_bsl;
            tmp3(fdr(Delay(i).F(b).TrType(tt).wplr_ltw_bsl_p)>0.05) = 0;
            imagesc(squeeze(tmp3));colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Baseline Corrected LTW');
            
            subplot(224)
            tmp4 = Delay(i).F(b).TrType(tt).wplr_ltw_bsl;
            tmp4(fdr(Delay(i).F(b).TrType(tt).wplr_ltw_bsl_p)>0.05) = 0;
            imagesc(squeeze(tmp4)); colorbar; caxis([-0.3,0.3]);
            set(gca, 'XTick',1:numel(rew_elecs), 'XTickLabel',rew_elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:numel(rew_elecs), 'YTickLabel',rew_elecs, 'YTickLabelRotation',45);
            title('Baseline Corrected MTW');
            flabel = strcat(num2str(Delay(i).length),'s',{' '},'Delay',{' '},freqs(b),...
                {' '},trTypes(tt),{' '},'Rew Regions WPLI Thresholded Connectivity');
            suptitle(flabel);
            keyboard;
            if(saveFigs)
                saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
            end
            %          
            %keyboard
            close;
        end
    end
end       

%% Combine data across all Delays for GLM input 
tic
for b = 1:6
    for tt = 1:2
        for tw = 1:3
            tot_ch_all = [];
            tot_ch_r = [];
            tot_wpl = [];
            tot_wplr = [];
            sesCnt = 0;
            sesCnt_r = 0;
            dcnt = 0;
            for i = 1:length(Delay)
                dcnt = dcnt+1;
                if(Delay(i).length==2) % exclude 2s delay
                    continue;
                end
                numS(dcnt) = size(Delay(i).F(b).TW(tw).TrType(tt).wpltmp,1);
                numS_r(dcnt) = size(Delay(i).F(b).TW(tw).TrType(tt).wplrtmp,1);
                sesCnt = sesCnt + numS(dcnt);
                sesCnt_r = sesCnt_r + numS_r(dcnt);
                tot_ch_all = vertcat(tot_ch_all,Delay(i).hr_ch_all');
                tot_ch_r = vertcat(tot_ch_r,Delay(i).hr_ch_allR');

                tot_wpl = cat(1,tot_wpl,Delay(i).F(b).TW(tw).TrType(tt).wpltmp);
                tot_wplr = cat(1,tot_wplr,Delay(i).F(b).TW(tw).TrType(tt).wplrtmp);
            end
            keyboard;
            all_ch_tot{b,tt,tw} = tot_ch_all;
            r_ch_tot{b,tt,tw} = tot_ch_r;
            wpl_tot{b,tt,tw} = tot_wpl;
            wplr_tot{b,tt,tw} = tot_wplr;

            d_1 = zeros(length(all_ch_tot{b,tt,tw}),1);
            d_2 = zeros(length(all_ch_tot{b,tt,tw}),1);
            d_3 = zeros(length(all_ch_tot{b,tt,tw}),1);
            d_4 = zeros(length(all_ch_tot{b,tt,tw}),1);
            sidx1 = numS(1);
            sidx2 = (numS(1)+numS(2));
            sidx3 = sidx2+numS(3);
            sidx4 = sidx3+numS(4);
            d_1(1:sidx1) = 1;
            d_2((sidx1+1):sidx2) = 1;
            d_3((sidx2+1):sidx3) = 1;
            d_4((sidx3+1):sidx4) = 1;
            glm_mat_wpl{b,tt,tw} = [d_1,d_2,d_3,d_4];
            
            d_1 = zeros(length(r_ch_tot{b,tt,tw}),1);
            d_2 = zeros(length(r_ch_tot{b,tt,tw}),1);
            d_3 = zeros(length(r_ch_tot{b,tt,tw}),1);
            d_4 = zeros(length(r_ch_tot{b,tt,tw}),1);
            sidx1 = numS_r(1);
            sidx2 = (numS_r(1)+numS_r(2));
            sidx3 = sidx2+numS_r(3);
            sidx4 = sidx3+numS_r(4);
            d_1(1:sidx1) = 1;
            d_2((sidx1+1):sidx2) = 1;
            d_3((sidx2+1):sidx3) = 1;
            d_4((sidx3+1):sidx4) = 1;
            glm_mat_wplr{b,tt,tw} = [d_1,d_2,d_3,d_4];
        end
    end
    
    % Concatenate HR-LR difference (outside of inner loop)
    tot_ch_r2 = [];
    tot_wplr2 = [];
    sesCnt_r2 = 0;
    dcnt = 0;
    for i = 1:length(Delay)
        dcnt = dcnt+1;
        if(Delay(i).length==2) % exclude 2s delay
            continue;
        end
        numS_r2(dcnt) = size(Delay(i).F(b).wplr_hrlrdf,1);
        sesCnt_r2 = sesCnt_r2 + numS_r2(dcnt);
        tot_ch_r2 = vertcat(tot_ch_r2,Delay(i).hr_ch_allR');
        tot_wplr2 = cat(1,tot_wplr2,Delay(i).F(b).wplr_hrlrdf);
    end
%     all_ch_tot2{b,tt,tw} = tot_ch_all2;
    r_ch_tot2{b} = tot_ch_r2;
    wplr_tot2{b} = tot_wplr2;

%     d_1 = zeros(length(all_ch_tot2{b,tt,tw}),1);
%     d_2 = zeros(length(all_ch_tot2{b,tt,tw}),1);
%     d_3 = zeros(length(all_ch_tot2{b,tt,tw}),1);
%     d_4 = zeros(length(all_ch_tot2{b,tt,tw}),1);
%     sidx1 = numS(1);
%     sidx2 = (numS(1)+numS(2));
%     sidx3 = sidx2+numS(3);
%     sidx4 = sidx3+numS(4);
%     d_1(1:sidx1) = 1;
%     d_2((sidx1+1):sidx2) = 1;
%     d_3((sidx2+1):sidx3) = 1;
%     d_4((sidx3+1):sidx4) = 1;
%     glm_mat_wpl_df{b,tt,tw} = [d_1,d_2,d_3,d_4];

    d_1b = zeros(length(r_ch_tot2{b}),1);
    d_2b = zeros(length(r_ch_tot2{b}),1);
    d_3b = zeros(length(r_ch_tot2{b}),1);
    d_4b = zeros(length(r_ch_tot2{b}),1);
    sidx1 = numS_r2(1);
    sidx2 = (numS_r2(1)+numS_r2(2));
    sidx3 = sidx2+numS_r2(3);
    sidx4 = sidx3+numS_r2(4);
    d_1b(1:sidx1) = 1;
    d_2b((sidx1+1):sidx2) = 1;
    d_3b((sidx2+1):sidx3) = 1;
    d_4b((sidx3+1):sidx4) = 1;
    glm_mat_wplr_df{b} = [d_1b,d_2b,d_3b,d_4b];
end
toc
disp('data combined across delays-GLM mat stored');            

%% Connectivity Behavioral Correlation (GLM)
% Correlate 12x12 reward reginal connectivity with %hr choice
for b = 1:6 
    for tt = 1:2
        for tw = 1:3
           for e =1:32 
               for e2 =1:32
                   if e == e2
                       continue;
                   end
                   L = length(wpl_tot{b,tt,tw}(:,e,e2));
                   glm_mat = [wpl_tot{b,tt,tw}(:,e,e2),glm_mat_wpl{b,tt,tw}(1:L,:)];
                   keyboard;
                   % with delay factor
                   [~,~,glm_stats{b,tt,tw,e,e2}] = glmfit(glm_mat, cell2mat(all_ch_tot{b,tt,tw}(1:L,:)));
                    
                   % without delay as factor
                   [~,~,glm_stats2{b,tt,tw,e,e2}] = glmfit(wpl_tot{b,tt,tw}(:,e,e2), cell2mat(all_ch_tot{b,tt,tw}(1:L,:)));

                   glm_beta{b,tt,tw,e,e2} = glm_stats{b,tt,tw,e,e2}.beta(2);

                   glm_p{b,tt,tw,e,e2} = glm_stats{b,tt,tw,e,e2}.p(2);

                   glm_beta2{b,tt,tw,e,e2} = glm_stats2{b,tt,tw,e,e2}.beta(2);

                   glm_p2{b,tt,tw,e,e2} = glm_stats2{b,tt,tw,e,e2}.p(2);
                end
           end
           for e = 1:12
               for e2 = 1:12
                   if e == e2
                       continue;
                   end
                   L = length(wplr_tot{b,tt,tw}(:,e,e2));
                   glm_mat_r = [wplr_tot{b,tt,tw}(:,e,e2),glm_mat_wplr{b,tt,tw}(1:L,:)];
                   [~,~,glm_stats_r{b,tt,tw,e,e2}] = glmfit(glm_mat_r, cell2mat(r_ch_tot{b,tt,tw}(1:L,:)));
                   [~,~,glm_stats2_r{b,tt,tw,e,e2}] = glmfit(wplr_tot{b,tt,tw}(:,e,e2), cell2mat(r_ch_tot{b,tt,tw}(1:L,:)));
                   
                   glm_beta_r{b,tt,tw,e,e2} = glm_stats_r{b,tt,tw,e,e2}.beta(2);
                   glm_p_r{b,tt,tw,e,e2} = glm_stats_r{b,tt,tw,e,e2}.p(2);
                   glm_beta2_r{b,tt,tw,e,e2} = glm_stats2_r{b,tt,tw,e,e2}.beta(2);
                   glm_p2_r{b,tt,tw,e,e2} = glm_stats2_r{b,tt,tw,e,e2}.p(2);
               end
           end
        end
    end
    
    % Calculate GLM  for HR-LR (2-1) difference across 12 rew regions
    for e = 1:12
       for e2 = 1:12
           if e == e2
               continue;
           end
           L = length(wplr_tot2{b}(:,e,e2));
           glm_mat_r_df = [wplr_tot2{b}(:,e,e2),glm_mat_wplr_df{b}(1:L,:)];
           [~,~,glm_stats_r_df{b,e,e2}] = glmfit(glm_mat_r_df, cell2mat(r_ch_tot2{b}(1:L,:)));
           [~,~,glm_stats2_r_df{b,e,e2}] = glmfit(wplr_tot2{b}(:,e,e2), cell2mat(r_ch_tot2{b}(1:L,:)));

           glm_beta_r_df{b,e,e2} = glm_stats_r_df{b,e,e2}.beta(2);
           glm_p_r_df{b,e,e2} = glm_stats_r_df{b,e,e2}.p(2);
           glm_beta2_r_df{b,e,e2} = glm_stats2_r_df{b,e,e2}.beta(2);
           glm_p2_r_df{b,e,e2} = glm_stats2_r_df{b,e,e2}.p(2);
       end
   end
end
disp('GLM calculated');
        
%% Plot GLM p and beta values with and without thresholding/FDR correction for 32 regions
evec = [28,27,26,25,31,32,29,30,6,5,1,20]; 
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];

elecs = ["NAcSh","A33","A24a","A24b","NAcCr","VMS","DMS","M1p","V1","V1","DorSub",...
    "DG","PPCx","CA1","CA3","STN","DLS","DLS","CE Am","BLA","A30c","A29c","MD Th",...
    "CE Th","vOFC","A32V","A32D","M2","L OFC","vAI","LFC","M1a"];
% plotting glm beta values 

for b = 4 
    for tw = 2    
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2 % plot only for middle time window 
            subplot(2,2,ttcnt)
            imagesc(squeeze(cell2mat(glm_beta(b,tt,tw,:,:)))); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:32, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:32, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            title(tlbl, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            colorbar;
            caxis([-2 2]);

            subplot(2,2,ttcnt2)
            imagesc(squeeze(cell2mat(glm_beta2(b,tt,tw,:,:)))); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).lr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).lr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:32, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:32, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl2 = strcat(TrType(tt),{' '},'No Delay Factored');
            title(tlbl2, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            colorbar;
            caxis([-2 2]);
            ttcnt=ttcnt+2;
            ttcnt2 = ttcnt2+2;
        end
        flabel = strcat(freqs(b),{' '},'GLM Connectivity Behavioral Correlation');
        suptitle(flabel);
    end
end

%% P value thresholded plots (FDR corrected) for 32 regions 

TrType = ["Low Reward","High Reward"];
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
elecs = ["NAcSh","A33","A24a","A24b","NAcCr","VMS","DMS","M1p","V1","V1","DorSub",...
    "DG","PPCx","CA1","CA3","STN","DLS","DLS","CE Am","BLA","A30c","A29c","MD Th",...
    "CE Th","vOFC","A32V","A32D","M2","L OFC","vAI","LFC","M1a"];

for b = 4 
    for tw = 2 % plot for MTW(0-1s RewL)
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2
            glm_p_c(b,tt,tw,:,:)= fdr(squeeze(cell2mat(glm_p(b,tt,tw,:,:))));
            glm_p2_c(b,tt,tw,:,:)= fdr(squeeze(cell2mat(glm_p2(b,tt,tw,:,:))));

            glm_beta_c=squeeze(cell2mat(glm_beta(b,tt,tw,:,:)));
            glm_beta_c(glm_p_c(b,tt,tw,:,:)>0.05)=0;

            glm_beta2_c=squeeze(cell2mat(glm_beta2(b,tt,tw,:,:)));
            glm_beta2_c(glm_p2_c(b,tt,tw,:,:)>0.05)=0;

            subplot(2,2,ttcnt)
            imagesc(glm_beta_c); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:32, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:32, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            title(tlbl, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            colorbar;
            caxis([-2 2]);

            subplot(2,2,ttcnt2)
            imagesc(glm_beta2_c); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).lr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).lr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:32, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:32, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl2 = strcat(TrType(tt),{' '},'No Delay Factored');
            title(tlbl2, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            colorbar;
            caxis([-2 2]);
            ttcnt=ttcnt+2;
            ttcnt2 = ttcnt2+2;
        end
       flabel = strcat(freqs(b),{' '},'FDR Corrected GLM Connectivity Behavioral Correlation');
       suptitle(flabel);
    end
end
%% Plot GLM p and beta values with and without thresholding/FDR correction for 12 rew regions
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
evec = [28,27,26,25,31,32,29,30,6,5,1,20]; 
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
% plotting glm beta values 
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
TrType = ["Low Reward","High Reward"];
saveFigs =1;
for b = [4 6] 
    for tw = 1    
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2 % plot only for middle time window 
            subplot(2,2,ttcnt)
            imagesc(squeeze(cell2mat(glm_beta_r(b,tt,tw,:,:)))); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            title(tlbl, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
%             colorbar;
            caxis([-2 2]);

            subplot(2,2,ttcnt2)
            imagesc(squeeze(cell2mat(glm_beta2_r(b,tt,tw,:,:)))); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).lr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).lr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl2 = strcat(TrType(tt),{' '},'No Delay Factored');
            title(tlbl2, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            colorbar;
            caxis([-2 2]);
            ttcnt=ttcnt+2;
            ttcnt2 = ttcnt2+2;
        end
        flabel = strcat(freqs(b),{' '},'GLM Connectivity Behavioral Correlation');
        suptitle(flabel);
        if(saveFigs)
            saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
end

%% P value thresholded plots (FDR corrected) for 12 rew regions 

TrType = ["Low Reward","High Reward"];
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];

for b = [4]
    for tw = 2 % plot for MTW(0-1s RewL)
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2
            tmp = cellfun(@isempty, glm_p_r);
            glm_p_r(tmp==1) = {NaN}; % set empty to 0
            tmp2 = cellfun(@isempty, glm_p2_r);
            glm_p2_r(tmp2==1) = {NaN}; % set empty to 0
            
            tmp3 = cellfun(@isempty, glm_beta_r);
            glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
            tmp4 = cellfun(@isempty, glm_beta2_r);
            glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0
            
            glm_p_cR(b,tt,tw,:,:)= fdr(squeeze(cell2mat(glm_p_r(b,tt,tw,:,:))));
            glm_p2_cR(b,tt,tw,:,:)= fdr(squeeze(cell2mat(glm_p2_r(b,tt,tw,:,:))));
        end


            glm_beta2_cR1=squeeze(cell2mat(glm_beta2_r(b,1,tw,:,:)));
            glm_beta2_cR1(glm_p2_cR(b,1,tw,:,:)>0.05)=0;
            glm_beta2_cR2=squeeze(cell2mat(glm_beta2_r(b,2,tw,:,:)));
            glm_beta2_cR2(glm_p2_cR(b,2,tw,:,:)>0.05)=0;

            
figure
            glm_beta2_cR1(isnan(glm_beta2_cR1)) = 0;
            imagesc(glm_beta2_cR1); % Display correlation matrix as an image
            set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl2 = strcat(TrType(1));
            title(tlbl2, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            caxis([-2 2]);
            colorbar
            axis square

        
%            subplot(1,3,2)
figure;
            glm_beta2_cR2(isnan(glm_beta2_cR2)) = 0;
            imagesc(glm_beta2_cR2); % Display correlation matrix as an image
            set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl2 = strcat(TrType(2));
            title(tlbl2, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            caxis([-2 2]);
                        colorbar;

            axis square

        
        
        
        flabel = strcat(freqs(b),{' '},'FDR Corrected GLM Connectivity Behavioral Correlation');
        suptitle(flabel);
        if(saveFigs)
           % saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
end

%%

TrType = ["Low Reward","High Reward"];
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];

for b = [4]
    for tw = 2 % plot for MTW(0-1s RewL)
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2
            tmp = cellfun(@isempty, glm_p_r);
            glm_p_r(tmp==1) = {NaN}; % set empty to 0
            tmp2 = cellfun(@isempty, glm_p2_r);
            glm_p2_r(tmp2==1) = {NaN}; % set empty to 0
            
            tmp3 = cellfun(@isempty, glm_beta_r);
            glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
            tmp4 = cellfun(@isempty, glm_beta2_r);
            glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0
            
            glm_p_cR(b,tt,tw,:,:)= squeeze(cell2mat(glm_p_r(b,tt,tw,:,:)));
            glm_p2_cR(b,tt,tw,:,:)= squeeze(cell2mat(glm_p2_r(b,tt,tw,:,:)));

            glm_beta_cR=squeeze(cell2mat(glm_beta_r(b,tt,tw,:,:)));
            glm_beta_cR(glm_p_cR(b,tt,tw,:,:)>0.05)=0;

            glm_beta2_cR=squeeze(cell2mat(glm_beta2_r(b,tt,tw,:,:)));
            glm_beta2_cR(glm_p2_cR(b,tt,tw,:,:)>0.05)=0;

            subplot(2,2,ttcnt)
            glm_beta_cR(isnan(glm_beta_cR)) = 0;
            imagesc(glm_beta_cR); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            title(tlbl, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
%             colorbar;
            caxis([-.06 .06]);

            subplot(2,2,ttcnt2)
            glm_beta2_cR(isnan(glm_beta2_cR)) = 0;
            imagesc(glm_beta2_cR); % Display correlation matrix as an image
            %     set(gca, 'XTick', 1:size(Delay(i).lr_P,1)); % center x-axis ticks on bins
            %     set(gca, 'YTick', 1:size(Delay(i).lr_P,2)); % center y-axis ticks on bins
            %     set(gca, 'XTickLabel', elecs); % set x-axis labels
            %     set(gca, 'YTickLabel', freqs); % set y-axis labels
            set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
            set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
            tlbl2 = strcat(TrType(tt),{' '},'No Delay Factored');
            title(tlbl2, 'FontSize', 10); % set title
            colormap('jet'); % Choose jet or any other color scheme
            colorbar;
            caxis([-.06 .06]);
            ttcnt=ttcnt+2;
            ttcnt2 = ttcnt2+2;
        end
        flabel = strcat(freqs(b),{' '},'FDR Corrected GLM Connectivity Behavioral Correlation');
        suptitle(flabel);
        if(saveFigs)
            saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
end



%% Plot GLM p and beta values with and without thresholding/FDR correction for 12 rew regions HR-LR Difference plots
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
evec = [28,27,26,25,31,32,29,30,6,5,1,20]; 
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
% plotting glm beta values 
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
TrType = ["Low Reward","High Reward"];
saveFigs =1;
for b = [4 6] 
    figure;
    subplot(1,2,1)
    imagesc(squeeze(cell2mat(glm_beta_r_df(b,:,:)))); % Display correlation matrix as an image
    %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
    %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
    %     set(gca, 'XTickLabel', elecs); % set x-axis labels
    %     set(gca, 'YTickLabel', freqs); % set y-axis labels
    set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
    set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
    tlbl = strcat('HR-LR Difference',{' '},'Delay Factored');
    title(tlbl, 'FontSize', 10); % set title
    colormap('jet'); % Choose jet or any other color scheme
%             colorbar;
    caxis([-.04 .04]);

    subplot(1,2,2)
    imagesc(squeeze(cell2mat(glm_beta2_r_df(b,:,:)))); % Display correlation matrix as an image
    %     set(gca, 'XTick', 1:size(Delay(i).lr_P,1)); % center x-axis ticks on bins
    %     set(gca, 'YTick', 1:size(Delay(i).lr_P,2)); % center y-axis ticks on bins
    %     set(gca, 'XTickLabel', elecs); % set x-axis labels
    %     set(gca, 'YTickLabel', freqs); % set y-axis labels
    set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
    set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
    tlbl2 = strcat('HR-LR Difference',{' '},'No Delay Factored');
    title(tlbl2, 'FontSize', 10); % set title
    colormap('jet'); % Choose jet or any other color scheme
    colorbar;
    caxis([-.04 .04]);
    flabel = strcat(freqs(b),{' '},'HRLR Diff GLM Connectivity Behavioral Correlation');
    flab2 = strcat(freqs(b),'GLM Connectivity Behavioral Correlation');
    suptitle(flab2);
    set(gcf,'Position',[751 489 760 384]);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
end

%%

%% P value thresholded plots (FDR corrected) for 12 rew regions HR-LR Difference plots 

TrType = ["Low Reward","High Reward"];
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];

for b = [4 6]
    figure;
    tmp = cellfun(@isempty, glm_p_r_df);
    glm_p_r_df(tmp==1) = {NaN}; % set empty to 0
    tmp2 = cellfun(@isempty, glm_p2_r_df);
    glm_p2_r_df(tmp2==1) = {NaN}; % set empty to 0

    tmp3 = cellfun(@isempty, glm_beta_r_df);
    glm_beta_r_df(tmp3==1) = {NaN}; % set empty to 0
    tmp4 = cellfun(@isempty, glm_beta2_r_df);
    glm_beta2_r_df(tmp4==1) = {NaN}; % set empty to 0

    glm_p_cR_df(b,:,:)= squeeze(cell2mat(glm_p_r_df(b,:,:)));
    glm_p2_cR_df(b,:,:)= squeeze(cell2mat(glm_p2_r_df(b,:,:)));

    glm_beta_cR_df=squeeze(cell2mat(glm_beta_r_df(b,:,:)));
    glm_beta_cR_df(fdr(glm_p_cR_df(b,:,:))>0.05)=0;

    glm_beta2_cR_df=squeeze(cell2mat(glm_beta2_r_df(b,:,:)));
    glm_beta2_cR_df(fdr(glm_p2_cR_df(b,:,:))>0.05)=0;

    subplot(1,2,1)
    glm_beta_cR_df(isnan(glm_beta_cR_df)) = 0;
    imagesc(glm_beta_cR_df); % Display correlation matrix as an image
    %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
    %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
    %     set(gca, 'XTickLabel', elecs); % set x-axis labels
    %     set(gca, 'YTickLabel', freqs); % set y-axis labels
    set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
    set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
    tlbl = strcat('HR-LR Difference',{' '},'Delay Factored');
    title(tlbl, 'FontSize', 10); % set title
    colormap('jet'); % Choose color scheme
%             colorbar;
    caxis([-.06 .06]);

    subplot(1,2,2)
    glm_beta2_cR_df(isnan(glm_beta2_cR_df)) = 0;
    imagesc(glm_beta2_cR_df); % Display correlation matrix as an image
    %     set(gca, 'XTick', 1:size(Delay(i).lr_P,1)); % center x-axis ticks on bins
    %     set(gca, 'YTick', 1:size(Delay(i).lr_P,2)); % center y-axis ticks on bins
    %     set(gca, 'XTickLabel', elecs); % set x-axis labels
    %     set(gca, 'YTickLabel', freqs); % set y-axis labels
    set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
    set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
    tlbl2 = strcat('HR-LR Difference',{' '},'No Delay Factored');
    title(tlbl2, 'FontSize', 10); % set title
    colormap('jet'); % Choose color scheme
    colorbar;
    caxis([-.06 .06]);
    flabel = strcat(freqs(b),{' '},'HRLR Diff FDR Corrected GLM Connectivity Behavioral Correlation');
    flab2 = strcat(freqs(b),'FDR Corrected GLM Connectivity Behavioral Correlation');
    suptitle(flab2);
    set(gcf,'Position',[751 489 760 384]);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
end



%% Connectivity Correlation Graph network plotting FDR corrected 12 elecs
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\GraphNetworkPlots\FDR Corrected\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
saveFigs = 1; 
graph_nodes = [{'M2'},{'A32D'},{'A32V'},{'vOFC'},{'LFC'},{'ALM'},{'LOFC'},{'Ains'},{'VMS'},{'NAcC'},{'NAcS'},{'BLA'}];
for b = [4 6]
    for tw = 2 % plot for MTW(0-1s RewL)
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2
            tmp = cellfun(@isempty, glm_p_r);
            glm_p_r(tmp==1) = {NaN}; % set empty to 0
            tmp2 = cellfun(@isempty, glm_p2_r);
            glm_p2_r(tmp2==1) = {NaN}; % set empty to 0
            
            tmp3 = cellfun(@isempty, glm_beta_r);
            glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
            tmp4 = cellfun(@isempty, glm_beta2_r);
            glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0
            
            glm_p_cR(b,tt,tw,:,:)= fdr(squeeze(cell2mat(glm_p_r(b,tt,tw,:,:))));
            glm_p2_cR(b,tt,tw,:,:)= fdr(squeeze(cell2mat(glm_p2_r(b,tt,tw,:,:))));

            glm_beta_cR=squeeze(cell2mat(glm_beta_r(b,tt,tw,:,:)));
            glm_beta_cR(glm_p_cR(b,tt,tw,:,:)>0.05)=0;

            glm_beta2_cR=squeeze(cell2mat(glm_beta2_r(b,tt,tw,:,:)));
            glm_beta2_cR(glm_p2_cR(b,tt,tw,:,:)>0.05)=0;
            
            
            glm_beta_cR(isnan(glm_beta_cR)) = 0;
            tmpP = glm_beta_cR;
            tmpN = glm_beta_cR;
            tmpP(tmpP<0) = 0;
            tmpN(tmpN>0) = 0;
            tmpSp = sparse(tmpP); % store sparse matrix
            tmpSn = sparse(tmpN); % store sparse matrix
            tmpGp = graph(tmpSp);
            tmpGn = graph(tmpSn);
            % Positive Connectivity correlation delay factored
            figure;
            g = plot(tmpGp,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network ');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            % Negative Connectivity correlation delay factored
            figure;
            g = plot(tmpGn,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            glm_beta2_cR(isnan(glm_beta2_cR)) = 0;
            tmpP = glm_beta2_cR;
            tmpN = glm_beta2_cR;
            tmpP(tmpP<0) = 0;
            tmpN(tmpN>0) = 0;
            tmpSp = sparse(tmpP); % store sparse matrix
            tmpSn = sparse(tmpN); % store sparse matrix
            tmpGp = graph(tmpSp);
            tmpGn = graph(tmpSn);
            % Positive Connectivity correlation no delay factored
            figure;
            g = plot(tmpGp,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'No Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            % Negative Connectivity correlation no delay factored
            figure;
            g = plot(tmpGn,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'No Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            close all;
        end
    end
end
cd(pdir)
 

%% No FDR correction graph network plotting 12 elecs
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\GraphNetworkPlots\NoFDR\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
saveFigs = 1; 
graph_nodes = [{'M2'},{'A32D'},{'A32V'},{'vOFC'},{'LFC'},{'ALM'},{'LOFC'},{'Ains'},{'VMS'},{'NAcC'},{'NAcS'},{'BLA'}];
for b = [4 6]
    for tw = 2 % plot for MTW(0-1s RewL)
        ttcnt = 1;
        ttcnt2 = 2;
        figure;
        for tt = 1:2
            tmp = cellfun(@isempty, glm_p_r);
            glm_p_r(tmp==1) = {NaN}; % set empty to 0
            tmp2 = cellfun(@isempty, glm_p2_r);
            glm_p2_r(tmp2==1) = {NaN}; % set empty to 0
            
            tmp3 = cellfun(@isempty, glm_beta_r);
            glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
            tmp4 = cellfun(@isempty, glm_beta2_r);
            glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0
            
            glm_p_cR(b,tt,tw,:,:)= squeeze(cell2mat(glm_p_r(b,tt,tw,:,:)));
            glm_p2_cR(b,tt,tw,:,:)= squeeze(cell2mat(glm_p2_r(b,tt,tw,:,:)));

            glm_beta_cR=squeeze(cell2mat(glm_beta_r(b,tt,tw,:,:)));
            glm_beta_cR(glm_p_cR(b,tt,tw,:,:)>0.05)=0;

            glm_beta2_cR=squeeze(cell2mat(glm_beta2_r(b,tt,tw,:,:)));
            glm_beta2_cR(glm_p2_cR(b,tt,tw,:,:)>0.05)=0;
            
            glm_beta_cR(isnan(glm_beta_cR)) = 0;
            tmpP = glm_beta_cR;
            tmpN = glm_beta_cR;
            tmpP(tmpP<0) = 0;
            tmpN(tmpN>0) = 0;
            tmpSp = sparse(tmpP); % store sparse matrix
            tmpSn = sparse(tmpN); % store sparse matrix
            tmpGp = graph(tmpSp);
            tmpGn = graph(tmpSn);
            % Positive Connectivity correlation delay factored
            figure;
            g = plot(tmpGp,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network ');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            % Negative Connectivity correlation delay factored
            figure;
            g = plot(tmpGn,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            glm_beta2_cR(isnan(glm_beta2_cR)) = 0;
            tmpP = glm_beta2_cR;
            tmpN = glm_beta2_cR;
            tmpP(tmpP<0) = 0;
            tmpN(tmpN>0) = 0;
            tmpSp = sparse(tmpP); % store sparse matrix
            tmpSn = sparse(tmpN); % store sparse matrix
            tmpGp = graph(tmpSp);
            tmpGn = graph(tmpSn);
            % Positive Connectivity correlation no delay factored
            figure;
            g = plot(tmpGp,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'No Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            % Negative Connectivity correlation no delay factored
            figure;
            g = plot(tmpGn,'Layout','force');
            labelnode(g,1:12,graph_nodes);
            tlbl = strcat(TrType(tt),{' '},'No Delay Factored');
            %title(tlbl, 'FontSize', 10); % set title
            flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
            title(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir,char(flabel)),'fig');
            end
            close all;
        end
    end
end                
cd(pdir)
                
            
%% Connectivity Correlation Graph network plotting FDR corrected 12 elecs HR-LR Difference 
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\GraphNetworkPlots\FDR Corrected\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
saveFigs = 1; 
graph_nodes = [{'M2'},{'A32D'},{'A32V'},{'vOFC'},{'LFC'},{'ALM'},{'LOFC'},{'Ains'},{'VMS'},{'NAcC'},{'NAcS'},{'BLA'}];
for b = [4 6]
    tmp = cellfun(@isempty, glm_p_r_df);
    glm_p_r(tmp==1) = {NaN}; % set empty to 0
    tmp2 = cellfun(@isempty, glm_p2_r_df);
    glm_p2_r(tmp2==1) = {NaN}; % set empty to 0

    tmp3 = cellfun(@isempty, glm_beta_r_df);
    glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
    tmp4 = cellfun(@isempty, glm_beta2_r_df);
    glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0

    glm_p_cR_df(b,:,:)= fdr(squeeze(cell2mat(glm_p_r_df(b,:,:))));
    glm_p2_cR_df(b,:,:)= fdr(squeeze(cell2mat(glm_p2_r_df(b,:,:))));

    glm_beta_cR_df=squeeze(cell2mat(glm_beta_r_df(b,:,:)));
    glm_beta_cR_df(glm_p_cR_df(b,:,:)>0.05)=0;

    glm_beta2_cR_df=squeeze(cell2mat(glm_beta2_r_df(b,:,:)));
    glm_beta2_cR_df(glm_p2_cR_df(b,:,:)>0.05)=0;

    glm_beta_cR_df(isnan(glm_beta_cR_df)) = 0;
    tmpP = glm_beta_cR_df;
    tmpN = glm_beta_cR_df;
    tmpP(tmpP<0) = 0;
    tmpN(tmpN>0) = 0;
    tmpSp = sparse(tmpP); % store sparse matrix
    tmpSn = sparse(tmpN); % store sparse matrix
    tmpGp = graph(tmpSp);
    tmpGn = graph(tmpSn);
    % Positive Connectivity correlation delay factored
    figure;
    g = plot(tmpGp,'Layout','force');
    labelnode(g,1:12,graph_nodes);
    tlbl = strcat('HR-LR Difference',{' '},'Delay Factored');
    %title(tlbl, 'FontSize', 10); % set title
    flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network ');
    title(flabel);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
    % Negative Connectivity correlation delay factored
    figure;
    g = plot(tmpGn,'Layout','force');
    labelnode(g,1:12,graph_nodes);
    tlbl = strcat('HR-LR Difference',{' '},'Delay Factored');
    %title(tlbl, 'FontSize', 10); % set title
    flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
    title(flabel);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
    glm_beta2_cR_df(isnan(glm_beta2_cR_df)) = 0;
    tmpP = glm_beta2_cR_df;
    tmpN = glm_beta2_cR_df;
    tmpP(tmpP<0) = 0;
    tmpN(tmpN>0) = 0;
    tmpSp = sparse(tmpP); % store sparse matrix
    tmpSn = sparse(tmpN); % store sparse matrix
    tmpGp = graph(tmpSp); % store Positive graph 
    tmpGn = graph(tmpSn); % store negative graph
    % Positive Connectivity correlation no delay factored
    figure;
    g = plot(tmpGp,'Layout','force');
    labelnode(g,1:12,graph_nodes);
    tlbl = strcat('HR-LR Difference',{' '},'No Delay Factored');
    %title(tlbl, 'FontSize', 10); % set title
    flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network');
    title(flabel);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
    % Negative Connectivity correlation no delay factored
    figure;
    g = plot(tmpGn,'Layout','force');
    labelnode(g,1:12,graph_nodes);
    tlbl = strcat('HR-LR Difference',{' '},'No Delay Factored');
    %title(tlbl, 'FontSize', 10); % set title
    flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
    title(flabel);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
    close all;
end
cd(pdir)

            
%% Connectivity Correlation Graph network plotting NO FDR corrected 12 elecs HR-LR Difference 
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\GraphNetworkPlots\NoFDR\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
saveFigs = 1; 
graph_nodes = [{'M2'},{'A32D'},{'A32V'},{'vOFC'},{'LFC'},{'ALM'},{'LOFC'},{'Ains'},{'VMS'},{'NAcC'},{'NAcS'},{'BLA'}];
for b = [4 6]
    tmp = cellfun(@isempty, glm_p_r_df);
    glm_p_r(tmp==1) = {NaN}; % set empty to 0
    tmp2 = cellfun(@isempty, glm_p2_r_df);
    glm_p2_r(tmp2==1) = {NaN}; % set empty to 0

    tmp3 = cellfun(@isempty, glm_beta_r_df);
    glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
    tmp4 = cellfun(@isempty, glm_beta2_r_df);
    glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0

    glm_p_cR_df(b,:,:)= squeeze(cell2mat(glm_p_r_df(b,:,:)));
    glm_p2_cR_df(b,:,:)= squeeze(cell2mat(glm_p2_r_df(b,:,:)));

    glm_beta_cR_df=squeeze(cell2mat(glm_beta_r_df(b,:,:)));
    glm_beta_cR_df(glm_p_cR_df(b,:,:)>0.05)=0;

    glm_beta2_cR_df=squeeze(cell2mat(glm_beta2_r_df(b,:,:)));
    glm_beta2_cR_df(glm_p2_cR_df(b,:,:)>0.05)=0;

    glm_beta_cR_df(isnan(glm_beta_cR_df)) = 0;
    tmpP = glm_beta_cR_df;
    tmpN = glm_beta_cR_df;
    tmpP(tmpP<0) = 0;
    tmpN(tmpN>0) = 0;
    if(~isempty(tmpP))
        tmpSp = sparse(tmpP); % store sparse matrix
        tmpGp = graph(tmpSp);
    end
    if(~isempty(tmpN))
        tmpSn = sparse(tmpN); % store sparse matrix
        tmpGn = graph(tmpSn);
    end
      
    if(exist('tmpGp','var'))
        % Positive Connectivity correlation delay factored
        figure;
        g = plot(tmpGp,'Layout','force');
        labelnode(g,1:12,graph_nodes);
        tlbl = strcat('HR-LR Difference',{' '},'Delay Factored');
        %title(tlbl, 'FontSize', 10); % set title
        flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network ');
        title(flabel);
        if(saveFigs)
            saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
    if(exist('tmpGn','var'))
        % Negative Connectivity correlation delay factored
        figure;
        g = plot(tmpGn,'Layout','force');
        labelnode(g,1:12,graph_nodes);
        tlbl = strcat('HR-LR Difference',{' '},'Delay Factored');
        %title(tlbl, 'FontSize', 10); % set title
        flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
        title(flabel);
        if(saveFigs)
            saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
    clear tmpGn tmpGp tmpSn tmpSp
    glm_beta2_cR_df(isnan(glm_beta2_cR_df)) = 0;
    tmpP = glm_beta2_cR_df;
    tmpN = glm_beta2_cR_df;
    tmpP(tmpP<0) = 0;
    tmpN(tmpN>0) = 0;
    if(~isempty(tmpP))
        tmpSp = sparse(tmpP); % store sparse matrix
        tmpGp = graph(tmpSp); % store Positive graph   
    end
    if(~isempty(tmpN))
        tmpSn = sparse(tmpN); % store sparse matrix
        tmpGn = graph(tmpSn); % store negative graph
    end
    if(exist('tmpGp','var'))
        % Positive Connectivity correlation no delay factored
        figure;
        g = plot(tmpGp,'Layout','force');
        labelnode(g,1:12,graph_nodes);
        tlbl = strcat('HR-LR Difference',{' '},'No Delay Factored');
        %title(tlbl, 'FontSize', 10); % set title
        flabel = strcat(freqs(b),{' '},tlbl,{' '},'Positive Correlation Network');
        title(flabel);
        if(saveFigs)
            saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
    if(exist('tmpGn','var'))
        % Negative Connectivity correlation no delay factored
        figure;
        g = plot(tmpGn,'Layout','force');
        labelnode(g,1:12,graph_nodes);
        tlbl = strcat('HR-LR Difference',{' '},'No Delay Factored');
        %title(tlbl, 'FontSize', 10); % set title
        flabel = strcat(freqs(b),{' '},tlbl,{' '},'Negative Correlation Network');
        title(flabel);
        if(saveFigs)
            saveas(gcf,fullfile(pdir,char(flabel)),'fig');
        end
    end
    pause;
    %close all;
end
cd(pdir)


          

%% Load PRL data for conjunction analysis
prldir = 'F:\PRL_Processed(storage)\physProcessed\Dynamic\respLV3\WT\Connectivity\BehCorrleation\';
load(fullfile(prldir,'PRL_ConnBehCorr_data.mat'));

%% Conjunction Analysis 
pdir = 'E:\DD_PhysProcessed\Connectivity\Results\WPLIPlots\CorrelationPlots\GraphNetworkPlots\NoFDR\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
% Make imagesc plot showing significance across LR(delay factored) negative map, 
% HR(no delay factor) 
TrType = ["Low Reward","High Reward"];
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
saveFigs=1;
for b = [4 6]
    % Win Stay Connnectivity Rev Count Correlation Matrices
    glm_pNR_cR = squeeze(cell2mat(glm_pR_nR(b,3,3,:,:)));
    glm_betaNR_cR = squeeze(cell2mat(glm_betaR_nR(b,3,3,:,:)));
    glm_betaNR_cR(glm_pNR_cR>0.05)=0;
    
    
    glm_pNR_cR_tc = squeeze(cell2mat(glm_pR_nR(b,1,3,:,:)));
    glm_betaNR_cR_tc = squeeze(cell2mat(glm_betaR_nR(b,1,3,:,:)));
    glm_betaNR_cR_tc(glm_pNR_cR_tc>0.05)=0;
    
    
    tconj = glm_betaNR_cR;
    
    tconj(glm_pNR_cR>0.05)=0;
%     tconj(glm_pNR_cR_tc>0.05)=0;


%     subplot(2,2,tt)
%     imagesc(glm_betaNR_cR); % Display correlation matrix as an image
%     set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
%     set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
%     tlbl = strcat(TrType(tt),{' '},'Reversal Count');
%     title(tlbl, 'FontSize', 10); % set title
%     colormap('jet'); % Choose jet or any other color scheme
%     if tt ==2 || tt == 4
%         colorbar;
%     end
%     caxis([-.04 .04]);
    
    tmp = cellfun(@isempty, glm_p_r);
    glm_p_r(tmp==1) = {NaN}; % set empty to 0
    tmp2 = cellfun(@isempty, glm_p2_r);
    glm_p2_r(tmp2==1) = {NaN}; % set empty to 0

    tmp3 = cellfun(@isempty, glm_beta_r);
    glm_beta_r(tmp3==1) = {NaN}; % set empty to 0
    tmp4 = cellfun(@isempty, glm_beta2_r);
    glm_beta2_r(tmp4==1) = {NaN}; % set empty to 0

    
    % Delay Factored Low Reward
    glm_p_cR = squeeze(cell2mat(glm_p_r(b,1,2,:,:)));
    glm_beta_cR = squeeze(cell2mat(glm_beta_r(b,1,2,:,:)));
    glm_beta_cR(glm_p_cR>0.05)=0;

    % No delay factored High Reward     
    glm_p2_cR = squeeze(cell2mat(glm_p2_r(b,2,2,:,:)));
    glm_beta2_cR = squeeze(cell2mat(glm_beta2_r(b,2,2,:,:)));
    glm_beta2_cR(glm_p2_cR>0.05)=0;
    
    % No delay factored Low Reward     
    glm_p2_cR_lr = squeeze(cell2mat(glm_p2_r(b,1,2,:,:)));
    glm_beta2_cR_lr = squeeze(cell2mat(glm_beta2_r(b,1,2,:,:)));
    glm_beta2_cR_lr(glm_p2_cR_lr>0.05)=0;
    
    % LR- Delay Factored DD tw=2
%     tconj(glm_p_cR>0.05) = 0;
    
    % LR- No Delay Factored DD tw=2  % threshold significant values between
    tconj(glm_p2_cR_lr>0.05) = 0;    % LR HR no delay factored
    
    % HR- No Delay Factored DD tw=2 
    tconj(glm_p2_cR>0.05) = 0;
    
    tconj(isnan(tconj)) = 0;
    figure;
    imagesc(tconj); % Display correlation matrix as an image
    %     set(gca, 'XTick', 1:size(Delay(i).hr_P,1)); % center x-axis ticks on bins
    %     set(gca, 'YTick', 1:size(Delay(i).hr_P,2)); % center y-axis ticks on bins
    %     set(gca, 'XTickLabel', elecs); % set x-axis labels
    %     set(gca, 'YTickLabel', freqs); % set y-axis labels
    set(gca, 'XTick',1:12, 'XTickLabel',elecs, 'XTickLabelRotation',45);
    set(gca, 'YTick',1:12, 'YTickLabel',elecs, 'YTickLabelRotation',30);
%     tlbl = strcat(TrType(tt),{' '},'Delay Factored');
%     title(tlbl, 'FontSize', 10); % set title
    colormap('jet'); % Choose jet or any other color scheme
%             colorbar;
    caxis([-.06 .06]);
    flabel = strcat(freqs(b),{' '},'Cojunctivity Analysis across LR-DF,HR-nDF,WS-nR');
    title(flabel);
    if(saveFigs)
        saveas(gcf,fullfile(pdir,char(flabel)),'fig');
    end
end
cd(pdir)










                
                
                 