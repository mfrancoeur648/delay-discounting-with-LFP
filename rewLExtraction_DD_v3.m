%%% RewL script
% Delay Discounting v1 binning script store entire time window across
% electrodes and fbands(bins the cannonical bands)
% Version 2: Properly normalize rewL for each trial type and save data 
% Version 3: Store entire frequency window for OFC spectograms
clearvars;
%% TC/nTC/TE/nTE
clear all; 
close all;

rdir='E:\dd_processed\trType_labeled\';
sdir ='E:\DD_PhysProcessed\rewL_processedDataV3\';
if~exist(sdir,'dir')
   mkdir(sdir);
end 
cd(rdir)
ddir = dir;
ddir(1:2)=[];
ancnt=1;
saveVars=1;
saveFigs=0;
plotFigs=0;
c1=0;
sescnt=0;
% switch based on which is running
tic
for i = 1:length(ddir) % Iterate through all delay lengths 
    disp(ddir(i).name)         
    cd([rdir '\' ddir(i).name]);
    adir = dir;
    adir(1:2) = []; 
    delayL(i) = str2num(extractBefore(ddir(i).name,'msH')); %#ok<ST2NM>
    for j = 1:length(adir) %Iterate through all sessions(across all anis)
        an = adir(j).name;
        if(startsWith(an,'R')&&contains(an,'_DD_'))                
            rat = extractBefore(an,'_DD_');
            if ancnt>1
                if(saveVars)
                    if(exist('D','var'))
                        savedir = strcat(sdir,num2str(delayL(i)),'msDelay\',rname{ancnt-1},'\');
                        if~exist(savedir,'dir')
                           mkdir(savedir);
                        end 
                        fname = strcat(rname{ancnt-1},'_rewL_fstored.mat');
                        save(fullfile(char(savedir),char(fname)),'D','-v7.3');
                        clear D;
                        disp(rname{ancnt-1});
                    end
                end
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
            if(exist('epochdata.mat','file'))                        
                sescnt = sescnt+1;
                D.An_phys{sescnt,1}=rat;
                D.An_phys{sescnt,2}=delayL(i);
                if(exist('extracted_data.mat','file'))
                    load('extracted_data.mat');
                    date = extractBetween(an,'_DD_','_');
%                      keyboard;
                    %load('epochdata2.mat');
                    D.behm= behm;
                    if(exist('DD_rew.mat','file'))
                        c2=0;
                        TF=load('DD_rew.mat');
                        if isfield(TF,'DD_rew')
                            DD_rew2 = TF.DD_rew.*conj(TF.DD_rew); %get power
                        end
                        bl=nanmean(nanmean(DD_rew2(:,:,25:35,:),3),4);
                        bl=repmat(bl,[1,1,200,size(DD_rew2,4)]);
                        %keyboard;
                        try
                            DD_rew2_blc=DD_rew2./bl; %base-line normalization

                            DD_rew2=10*log10(DD_rew2); % log transforming the power.
                            DD_rew2_blc=10*log10(DD_rew2_blc); % log transforming the power.

                            bl=nanmean(nanmean(DD_rew2(:,:,25:35,:),3),4);
                            bl=repmat(bl,[1,1,200,size(DD_rew2,4)]);
                            DD_rew2=DD_rew2-bl; % zero-centering the data.

                            bl_blc=nanmean(nanmean(DD_rew2_blc(:,:,25:35,:),3),4);
                            bl_blc=repmat(bl_blc,[1,1,200,size(DD_rew2_blc,4)]);
                            DD_rew2_blc=DD_rew2_blc-bl_blc; % zero-centering the data.
                        catch ME 
                             if (strcmp(ME.identifier,'MATLAB:Array dimensions must match for binary array op.'))
                                 sescnt = sescnt-1;
                                 continue;
                             end
                        end
                        E = size(DD_rew2_blc,1);
                        L = size(DD_rew2_blc,4);
                        lb = [0,4,8,13,30,70];
                        hb = [4,8,13,30,50,120];
                        for e = 1:E % iterate through all electrodes
                            lrcnt = 0;
                            hrcnt = 0;
                            for tr = 1:L % iterate through all trials within session
                                % take mean across all 6 cannonical freq windows 
                                rewL_tr{e,tr} = squeeze(mean(DD_rew2_blc(e,:,:,tr),2,'omitnan'));
                                % Store low and high reward vectors
                                % across all trials 
                                np2v(tr) = isempty(behm.res_2{tr});
                                np4v(tr) = isempty(behm.res_4{tr});
                                if(np2v(tr)==1)
                                    lrcnt = lrcnt+1;
                                    rewL_lr(lrcnt,:,:) = rewL_tr{e,tr};
                                elseif(np4v(tr)==1)
                                    hrcnt = hrcnt+1;
                                    rewL_hr(hrcnt,:,:) = rewL_tr{e,tr};
                                end
                            end
                            % Store 3 dim data for high and low reward
                            % average across trials 
                            nanvec = NaN(1,200);
                            if exist('rewL_lr','var')
                                D.rewL_lr{sescnt,e} = squeeze(mean(rewL_lr,1,'omitnan')); %#ok<*UDIM>
                            else
                                D.rewL_lr{sescnt,e} = nanvec;
                            end
                            if exist('rewL_hr','var')
                                D.rewL_hr{sescnt,e} = squeeze(mean(rewL_hr,1,'omitnan'));
                            else 
                                D.rewL_hr{sescnt,e} = nanvec;
                            end
                            clear rewL_hr rewL_lr
                        end
                    end
                end
            end
        end
    end
end
toc
disp('DD version 2 rewL data extracted'); 

                        
                        
                        
      
                                    
                                    
                                    