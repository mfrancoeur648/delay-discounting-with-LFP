% Delay Discounting Mean Power plotting rewLplottingDD_v1
% Version 1 - rewL Imagesc plots yaxis sessions
% for electrode 29( to check if data appears across all sessions
% Version 2 - rewL imagesc plots y xis electrodes x axis time- mean power
% across sessions animals and for hr power, lr power
% Version 3 - Exp plotting script to compare hr lr and difference
% Version 4 zscEdit- rewL imagesc plots - zscored for rew regions
% Version 5 zscEdit - shaded error bar plots for OFC region LR vs HR across
% low medium and high delay ( paper plots grouped)
% Version 7: Rebin data- Average across all sessions pooled across all
% animals 
% Version 9: Paper specific plotting with DD zscored data 
clearvars;
close all;
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
 % Early Time Window(-300-0ms)
etw1 = -300;
etw2 = 0;

etwa1 = etw1-30;
etwa2 = etw2+30;
tt_e = t(t>etwa1 & t<etwa2);
etvec = find(t==min(tt_e)):find(t==max(tt_e));
ebinidx1 = min(etvec);
ebinidx2 = max(etvec);

% Late Time Window(0-1000ms)
ltw1 = 0000;
ltw2 = 1000;

ltwa1 = ltw1-30;
ltwa2 = ltw2+30;
tt_l = t(t>ltwa1 & t<ltwa2);
ltvec = find(t==min(tt_l)):find(t==max(tt_l));
lbinidx1 = min(ltvec);
lbinidx2 = max(ltvec);

% Base Line Time Window(0-1000ms)
bstw1 = -2500;
bstw2 = -1500;

bstwa1 = bstw1-30;
bstwa2 = bstw2+30;
tt_b = t(t>bstwa1 & t<bstwa2);
bslvec = find(t==min(tt_b)):find(t==max(tt_b));
blbinidx1 = min(bslvec);
blbinidx2 = max(bslvec);

% Base Line Time Window 2 (0-1000ms)
bstw1 = -2500;
bstw2 = -2000;

bstwa1 = bstw1-30;
bstwa2 = bstw2+30;
tt_b = t(t>bstwa1 & t<bstwa2);
bslvec2 = find(t==min(tt_b)):find(t==max(tt_b));
blbinidx1 = min(bslvec);
blbinidx2 = max(bslvec);



saveVars=1;
saveFigs=1;
plotFigs=0;
c1=0;
sescnt=0;
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
% switch based on which is running
tic
for i = 1:length(ddir) % Iterate through all delay lengths
    disp(ddir(i).name)
    cd([rdir '\' ddir(i).name]);
    adir = dir;
    adir(1:2) = [];
    ancnt=0;
    scnt = 0;
    Delay(i).length = str2num(extractBefore(ddir(i).name,'ms'));
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
            for e = 1:E
                for b = 1:6
                    numSes = size(D.rewL_hr_blc,1);
                    if numSes>1
                        for s = 1:numSes
                            nanvec = NaN(1,200);
                            if(length(cell2mat(D.rewL_hr_blc(s,e,b)))~=200)
                                D.rewL_hr_blc{s,e,b} = nanvec;
                            end
                            if(length(cell2mat(D.rewL_lr_blc(s,e,b)))~=200)
                                D.rewL_lr_blc{s,e,b} = nanvec;
                            end
%                             hrv(s,:) = zscore(D.rewL_hr_blc{s,e,b});
%                            lrv(s,:) = zscore(D.rewL_lr_blc{s,e,b});
                   %         dfv(s,:) = hrv(s,:) - lrv(s,:);
                            
                            if(Delay(i).length==2)
                                hrv(s,:) = (D.rewL_hr_blc{s,e,b}-...
                                    mean(D.rewL_hr_blc{s,e,b}(bslvec2)))./std(D.rewL_hr_blc{s,e,b}(bslvec2));
                                lrv(s,:) = (D.rewL_lr_blc{s,e,b}-...
                                    mean(D.rewL_lr_blc{s,e,b}(bslvec2)))./std(D.rewL_lr_blc{s,e,b}(bslvec2));
                                dfv(s,:) = hrv(s,:) - lrv(s,:);
%                                 hrv(s,:) = zscore(D.rewL_hr_blc{s,e,b});
%                                 lrv(s,:) = zscore(D.rewL_lr_blc{s,e,b});
%                                 dfv(s,:) = hrv(s,:) - lrv(s,:);
                            else
                                hrv(s,:) = (D.rewL_hr_blc{s,e,b}-...
                                    mean(D.rewL_hr_blc{s,e,b}(bslvec)))./std(D.rewL_hr_blc{s,e,b}(bslvec));
                                lrv(s,:) = (D.rewL_lr_blc{s,e,b}-...
                                    mean(D.rewL_lr_blc{s,e,b}(bslvec)))./std(D.rewL_lr_blc{s,e,b}(bslvec));
                                dfv(s,:) = hrv(s,:) - lrv(s,:);
                            end
                            
                            Delay(i).rewL_hr{ancnt,s,e,b,:} = hrv(s,:);
                            Delay(i).rewL_lr{ancnt,s,e,b,:} =  lrv(s,:);
                            Delay(i).rewL_df{ancnt,s,e,b,:} = dfv(s,:);
                        end
                    elseif numSes==1
                        s=1;
                        nanvec = NaN(1,200);
                        if(length(cell2mat(D.rewL_hr_blc(s,e,b)))~=200)
                            D.rewL_hr_blc{s,e,b} = nanvec;
                        end
                        if(length(cell2mat(D.rewL_lr_blc(s,e,b)))~=200)
                            D.rewL_lr_blc{s,e,b} = nanvec;
                        end
                         hrv(s,:) = zscore(D.rewL_hr_blc{s,e,b});
                         lrv(s,:) = zscore(D.rewL_lr_blc{s,e,b});
                         hrv(s,:) = (D.rewL_hr_blc{s,e,b}-mean(D.rewL_hr_blc{s,e,b}(bslvec)))./std(D.rewL_hr_blc{s,e,b}(bslvec));
                         lrv(s,:) = (D.rewL_lr_blc{s,e,b}-mean(D.rewL_lr_blc{s,e,b}(bslvec)))./std(D.rewL_lr_blc{s,e,b}(bslvec));
                            
                         dfv(s,:) = hrv(s,:) - lrv(s,:);
                        
                         Delay(i).rewL_hr{ancnt,s,e,b,:} = hrv(s,:);
                         Delay(i).rewL_lr{ancnt,s,e,b,:} = lrv(s,:);
                         Delay(i).rewL_df{ancnt,s,e,b,:} = dfv(s,:);
                    end
                end
            end
        end
    end
end
toc
disp('DD rewL plotting preprocessed');

%% Concatenate data across sessions
tic
for i = 1:length(Delay)-1 % iterate delays skip 40s delay(not enough hr tr)
    numAnis = size(Delay(i).rewL_hr,1);
    numSess = size(Delay(i).rewL_hr,2);
    numE = size(Delay(i).rewL_hr,3);
    numB = size(Delay(i).rewL_hr,4);
    for e = 1:numE %iterate elecs
        for b = 1:numB %iterate fbands
            % allocate for each delay (new concat for each elec/ freq)
            hrvec = [];
            lrvec = [];
            dfvec = [];
            for a = 1:numAnis  %iterate over animals
                for s = 1:numSess  %iterate sessions
                    % store data across all animals and sessions
                    hrvec = vertcat(hrvec,Delay(i).rewL_hr{a,s,e,b});
                    lrvec = vertcat(lrvec,Delay(i).rewL_lr{a,s,e,b});
                    dfvec = vertcat(dfvec,Delay(i).rewL_df{a,s,e,b});
                end
            end
         %   hr_vec{e,b} = hrvec;
          %  lr_vec{e,b} = lrvec;
           
          %df_vec{e,b} = dfvec;    
            
            Delay(i).hrvec{e,b} = hrvec;
            Delay(i).lrvec{e,b} = lrvec;
            Delay(i).dfvec{e,b} = dfvec;
            % Store mean value across all sessions
            Delay(i).rewL_hr_mn{e,b} = mean(Delay(i).hrvec{e,b},1,'omitnan');
            Delay(i).rewL_lr_mn{e,b} = mean(Delay(i).lrvec{e,b},1,'omitnan');
            Delay(i).rewL_df_mn{e,b} = mean(Delay(i).dfvec{e,b},1,'omitnan');

            Delay(i).rewL_hr_sem{e,b,:} = std(Delay(i).hrvec{e,b},1,'omitnan')/...
                sqrt(size(Delay(i).hrvec{e,b},1));
            Delay(i).rewL_lr_sem{e,b,:} = std(Delay(i).lrvec{e,b},1,'omitnan')/...
                sqrt(size(Delay(i).lrvec{e,b},1));
            Delay(i).rewL_df_sem{e,b,:} = std(Delay(i).dfvec{e,b},1,'omitnan')/...
                sqrt(size(Delay(i).dfvec{e,b},1));

            tmp1 =  mean(Delay(i).hrvec{e,b}(:,lbinidx1:lbinidx2),2,'omitnan');
            tmp2 =  mean(Delay(i).lrvec{e,b}(:,lbinidx1:lbinidx2),2,'omitnan');
            tmp3 =  mean(Delay(i).dfvec{e,b}(:,lbinidx1:lbinidx2),2,'omitnan');
             
            % Store mean across Late Time Window (0-1000ms)
            Delay(i).hr_ltw_mn(e,b) = mean(tmp1,'omitnan');
            Delay(i).lr_ltw_mn(e,b) = mean(tmp2,'omitnan');
            Delay(i).df_ltw_mn(e,b) = mean(tmp3,'omitnan');
            
            Delay(i).hr_ltw_sem(e,b) = std(tmp1,'omitnan')/sqrt(length(tmp1));
            Delay(i).lr_ltw_sem(e,b) = std(tmp2,'omitnan')/sqrt(length(tmp2));
            Delay(i).df_ltw_sem(e,b) = std(tmp3,'omitnan')/sqrt(length(tmp3));
        end
    end
end
toc
disp('Data Concatenated across sessions');


%% Grouped plotting hr vs lr
close all;
pdir = 'E:\DD_PhysProcessed\rewL plots\spectAllEs\rewregions_zscored\lmh_grouped\';
if ~exist(pdir,'dir')
    mkdir(pdir)
end
pdir2 = 'E:\DD_PhysProcessed\rewL plots\paperplots\lmh_grouped_comps\updated\';
if ~exist(pdir2,'dir')
    mkdir(pdir2)
end

load('E:\DD_PhysProcessed\times_freqs.mat');
% Desired Time Window for plotting 

tw1 = -500;
tw2 = 2500;

twa1 = tw1-30;
twa2 = tw2+30;
tt = t(t>twa1 & t<twa2);
ttvec = find(t==min(tt)):find(t==max(tt));

% Load time and frequency vectors for plotting
tic
target_e = 29;

plotFigs=0;
plotFigs1 =0; % updated figure
plotFigs2 =1; % updated to include bars across all delay conditions
evec = [1,5,6,7,25,26,27,28,29];
freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
%elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
elecs = ["NAcSh","A33","A24a","A24b","NAcCr","VMS","DMS","M1p","V1","V1","DorSub",...
    "DG","PPCx","CA1","CA3","STN","DLS","DLS","CE Am","BLA","A30c","A29c","MD Th",...
    "CE Th","vOFC","A32V","A32D","M2","L OFC","vAI","LFC","M1a"];
twlab1 = tw1:((tw2-tw1+1)/length(ttvec)):tw2;
evec = [28,27,26,25,31,32,29,30,6,5,1,20]; 
%%

% 12 target rew regions 
% plot low med and high delay grouping imagesc plots
% for b = 1:6
%     ecnt=0;
%     for e = evec 
%         hr_all=[];
%         lr_all=[];
%         hr_ld=[];
%         lr_ld=[];
%         hr_md=[];
%         lr_md=[];
%         hr_hd=[];
%         lr_hd=[];
%         ecnt = ecnt+1;
%         for i = 1:length(Delay)-1  % Iterate through all delay lengths 
%             tmp_h = cell2mat(Delay(i).rewL_hr_mn(e,b)); %
%             tmp_l = cell2mat(Delay(i).rewL_lr_mn(e,b)); %
%             tmp_d = cell2mat(Delay(i).rewL_df_mn(e,b)); %
% 
%             tmp_h_sem = cell2mat(Delay(i).rewL_hr_sem(e,b)); %
%             tmp_l_sem = cell2mat(Delay(i).rewL_lr_sem(e,b)); %
%             tmp_d_sem = cell2mat(Delay(i).rewL_df_sem(e,b)); %
% 
%             D(i).length = Delay(i).length;
%             % store for each delay length
%             hr_zsc = tmp_h;
%             lr_zsc = tmp_l;
%             df_zsc = tmp_l;
%             
%             
%             % Remove bad animals
%             bad1 = isnan(hr_zsc(:,1));
%             hr_zsc(bad1,:) = [];
%             bad2 = isnan(lr_zsc(:,1));
%             lr_zsc(bad2,:) = [];
%             bad3 = isnan(df_zsc(:,1));
%             df_zsc(bad3,:) = [];
%             
%             
%             elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
%             evec = [28,27,26,25,31,32,29,30,6,5,1,20];
%             
%             % Take mean across Late Time Windows
%             hr_ltw_mn = Delay(i).hr_ltw_mn(evec,b); % mean(Delay(i).hr_ltw(:,ebinidx1:ebinidx2),2,'omitnan');
%             lr_ltw_mn = Delay(i).lr_ltw_mn(evec,b); % mean(lr_zsc(:,ebinidx1:ebinidx2),2,'omitnan');
%             df_ltw_mn = Delay(i).df_ltw_mn(evec,b); % mean(df_zsc(:,ebinidx1:ebinidx2),2,'omitnan');
% 
%             hr_ltw_sem = Delay(i).hr_ltw_sem(evec,b); % mean(Delay(i).hr_ltw(:,ebinidx1:ebinidx2),2,'omitnan');
%             lr_ltw_sem = Delay(i).lr_ltw_sem(evec,b); % mean(lr_zsc(:,ebinidx1:ebinidx2),2,'omitnan');
%             df_ltw_sem = Delay(i).df_ltw_sem(evec,b); % mean(df_zsc(:,ebinidx1:ebinidx2),2,'omitnan');            
%         end
% 
% 
% 
%         evec = [28,27,26,25,31,32,29,30,6,5,1,20];
%         
%         
% 
%     end
%%
bvec = [1,2,3,4,5,6];
evec = [28,27,26,25,31,32,29,30,6,5,1,20];
for b = 3
    if(plotFigs2) % only plot for beta
        for elec = evec
            %target_e=29;
    %        elec = find(evec == target_e); % plot for single chosen electrode
            figure;
            set(gcf,'Position', [-35   591   792   376]);
    %        twlab1 = -500:(3001/length(ttvec)):2500;
            twlab1=t;
            ttvec=1:200;
            colvec = ['b','r','g','m'];
            dvec = [1,3,5,6,2,4]; % Order Delays are stored in structure

            axl=[-2500 2500 -2 10]
            % shaded err bar plots across 3 delay conditions
            subplot(331)
            hold on;
            % 1ms shaded err bar plot
            shadedErrorBar(twlab1,Delay(dvec(1)).rewL_hr_mn{elec,b}(ttvec),Delay(dvec(1)).rewL_hr_sem{elec,b}(ttvec),...
                'lineprops',colvec(1),'patchSaturation',0.33);
            shadedErrorBar(twlab1,Delay(dvec(1)).rewL_lr_mn{elec,b}(ttvec),Delay(dvec(1)).rewL_lr_sem{elec,b}(ttvec),...
                'lineprops',colvec(2),'patchSaturation',0.33); 
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
            %xlabel('Time(ms)');
            ylabel('Zscored Power'); 
            title('500ms Delay')
                    axis(axl);

            %legend('HR','LR','Location','best'); 
            subplot(332)
            hold on;
            % 5ms shaded err bar plot
            shadedErrorBar(twlab1,Delay(dvec(2)).rewL_hr_mn{elec,b}(ttvec),Delay(dvec(2)).rewL_hr_sem{elec,b}(ttvec),...
                'lineprops',colvec(1),'patchSaturation',0.33);
            shadedErrorBar(twlab1,Delay(dvec(2)).rewL_lr_mn{elec,b}(ttvec),Delay(dvec(2)).rewL_lr_sem{elec,b}(ttvec),...
                'lineprops',colvec(2),'patchSaturation',0.33); 
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
            %xlabel('Time(ms)');
           % ylabel('Zscored Power'); 
            title('1S Delay')
                    axis(axl);
            %legend('HR','LR','Location','best');
            subplot(333)
            hold on;
            % 20 ms shaded err bar plot
            shadedErrorBar(twlab1,Delay(dvec(3)).rewL_hr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_hr_sem{elec,b}(ttvec),...
                'lineprops',colvec(1),'patchSaturation',0.33);
            shadedErrorBar(twlab1,Delay(dvec(3)).rewL_lr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_lr_sem{elec,b}(ttvec),...
                'lineprops',colvec(2),'patchSaturation',0.33); 
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
    %        xlabel('Time(ms)');
            %ylabel('Zscored Power'); 
            title('2S Delay')
      %      legend('HR','LR','Location','best');
                    axis(axl);

            subplot(334)
            hold on;
            % 20 ms shaded err bar plot
            shadedErrorBar(twlab1,Delay(dvec(4)).rewL_hr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_hr_sem{elec,b}(ttvec),...
                'lineprops',colvec(1),'patchSaturation',0.33);
            shadedErrorBar(twlab1,Delay(dvec(4)).rewL_lr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_lr_sem{elec,b}(ttvec),...
                'lineprops',colvec(2),'patchSaturation',0.33); 
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
    %        xlabel('Time(ms)');
            %ylabel('Zscored Power'); 
            title('5S Delay')
     %       legend('HR','LR','Location','best');
                    axis(axl);

            subplot(335)
            hold on;
            % 20 ms shaded err bar plot
            shadedErrorBar(twlab1,Delay(dvec(5)).rewL_hr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_hr_sem{elec,b}(ttvec),...
                'lineprops',colvec(1),'patchSaturation',0.33);
            shadedErrorBar(twlab1,Delay(dvec(5)).rewL_lr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_lr_sem{elec,b}(ttvec),...
                'lineprops',colvec(2),'patchSaturation',0.33); 
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
    %        xlabel('Time(ms)');
            %ylabel('Zscored Power'); 
            title('10S Delay')
    %        legend('HR','LR','Location','best');
                    axis(axl);

            subplot(336)
            hold on;
            % 20 ms shaded err bar plot
            shadedErrorBar(twlab1,Delay(dvec(6)).rewL_hr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_hr_sem{elec,b}(ttvec),...
                'lineprops',colvec(1),'patchSaturation',0.33);
            shadedErrorBar(twlab1,Delay(dvec(6)).rewL_lr_mn{elec,b}(ttvec),Delay(dvec(6)).rewL_lr_sem{elec,b}(ttvec),...
                'lineprops',colvec(2),'patchSaturation',0.33); 
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
    %        xlabel('Time(ms)');
            %ylabel('Zscored Power'); 
            title('20S Delay')
            legend('HR','LR','Location','best');
                    axis(axl);


            % traces plotted for each delay condition(l,m,h) LR, HR, Difference
            % plot
            colvec = ['r','r','c','c','k','k','y'];
            dvec = [1,3,5,6,2,4]; % Order Delays are stored in structure
            subplot(337) % LR plot
            tt=t;
            for d = [1,3,5]
                plot(tt,Delay(dvec(d)).rewL_lr_mn{elec,b}(ttvec),'Color',colvec(d)); hold on;
                plot(tt,Delay(dvec(d+1)).rewL_lr_mn{elec,b}(ttvec),'--','Color',colvec(d+1)); hold on;

            end
            title('Low Reward');
           % xlabel('Time(ms)');
            ylabel('Zscored Power');
                    axis(axl);
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
            plot([tt(1);tt(end)],[0;0],'k--');
            box off;

            %legend('LD','MD','HD');

            subplot(338) % HR plot
            for d = [1,3,5]
                plot(tt,Delay(dvec(d)).rewL_hr_mn{elec,b}(ttvec),'Color',colvec(d)); hold on;
               plot(tt,Delay(dvec(d+1)).rewL_hr_mn{elec,b}(ttvec),'--','Color',colvec(d+1)); hold on;

            end
            title('High Reward');
            xlabel('Time(ms)');
            %ylabel('Zscored Power');
                    axis(axl);
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
            plot([tt(1);tt(end)],[0;0],'k--');
            box off;                
            %legend('LD','MD','HD');

            subplot(339) % HR-LR difference plot
            for d = [1,3,5]
                plot(tt,Delay(dvec(d)).rewL_df_mn{elec,b}(ttvec),'Color',colvec(d)); hold on;
                plot(tt,Delay(dvec(d+1)).rewL_df_mn{elec,b}(ttvec),'--','Color',colvec(d+1)); hold on;

            end
            title('High Rew - Low Rew');
           % xlabel('Time(ms)');
            %ylabel('Zscored Power');
                    axis(axl);
            yl = ylim;
            plot([0,0,yl(1),yl(2)],'--k')
            plot([tt(1);tt(end)],[0;0],'k--');
            box off;                
            legend('.5s','1','2','5','10','20','Location','best');

            flabel = freqs(b);
            suptitle(flabel);
            %%%%%%

            % bar plots for individual electrode ETW and LTW across all delays
            % ETW plots - LR 
            %delayList= ["Low","Med","High"];
            delayList= ["0.5","1","2","5","10","20"];
    %
            figure;
            set (gcf,'Position',[680   302   250   676]);
            % LTW plots - LR 
            subplot(3,1,1)
            bars = [Delay(dvec(1)).hr_ltw_mn(elec,b),Delay(dvec(2)).hr_ltw_mn(elec,b),Delay(dvec(3)).hr_ltw_mn(elec,b),...
                Delay(dvec(4)).hr_ltw_mn(elec,b),Delay(dvec(5)).hr_ltw_mn(elec,b),Delay(dvec(6)).hr_ltw_mn(elec,b)];
            errs = [Delay(dvec(1)).hr_ltw_sem(elec,b),Delay(dvec(2)).hr_ltw_sem(elec,b),Delay(dvec(3)).hr_ltw_sem(elec,b),...
                Delay(dvec(4)).hr_ltw_sem(elec,b),Delay(dvec(5)).hr_ltw_sem(elec,b),Delay(dvec(6)).hr_ltw_sem(elec,b)];

            b1 = bar(bars,'grouped'); hold on;
            errorbar(bars, errs, 'r' , 'linestyle', 'none');

            xticklabels(delayList);
             xtickangle(90);
           % xlabel('Delay Lengths (S)');
            %ylabel('Zscored Power');
            title('High Reward');
            hold off; box off


            subplot(3,1,2) % HR ltw plot
            bars = [Delay(dvec(1)).lr_ltw_mn(elec,b),Delay(dvec(2)).lr_ltw_mn(elec,b),Delay(dvec(3)).lr_ltw_mn(elec,b),...
                Delay(dvec(4)).lr_ltw_mn(elec,b),Delay(dvec(5)).lr_ltw_mn(elec,b),Delay(dvec(6)).lr_ltw_mn(elec,b)];
            errs = [Delay(dvec(1)).lr_ltw_sem(elec,b),Delay(dvec(2)).lr_ltw_sem(elec,b),Delay(dvec(3)).lr_ltw_sem(elec,b),...
                Delay(dvec(4)).lr_ltw_sem(elec,b),Delay(dvec(5)).lr_ltw_sem(elec,b),Delay(dvec(6)).lr_ltw_sem(elec,b)];
            b1 = bar(bars,'grouped'); hold on;
            errorbar(bars, errs, 'r' , 'linestyle', 'none');

            xticklabels(delayList);
             xtickangle(90);
            %xlabel('Delay Lengths');
            ylabel('Zscored Power');
            title('Low Reward');
            hold off; box off



            subplot(3,1,3)
           bars = [Delay(dvec(1)).df_ltw_mn(elec,b),Delay(dvec(2)).df_ltw_mn(elec,b),Delay(dvec(3)).df_ltw_mn(elec,b),...
                Delay(dvec(4)).df_ltw_mn(elec,b),Delay(dvec(5)).df_ltw_mn(elec,b),Delay(dvec(6)).df_ltw_mn(elec,b)];
            errs = [Delay(dvec(1)).df_ltw_sem(elec,b),Delay(dvec(2)).df_ltw_sem(elec,b),Delay(dvec(3)).df_ltw_sem(elec,b),...
                Delay(dvec(4)).df_ltw_sem(elec,b),Delay(dvec(5)).df_ltw_sem(elec,b),Delay(dvec(6)).df_ltw_sem(elec,b)];

            b1 = bar(bars,'grouped'); hold on;
            errorbar(bars, errs, 'r' , 'linestyle', 'none');
            xticklabels(delayList);
             xtickangle(90);
            xlabel('Delay Lengths (S)');
            %ylabel('Zscored Power');
            title('High Rew - Low Rew');
            hold off; box off
            flabel=strcat(freqs(b),{' '},'Electrode',num2str(elec),{' '},'Delay Comparisons');
            suptitle(flabel);
            if(saveFigs)
                saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
            end
            pause;
        end
    end
end
        
%% Bar plots across all rew related regions late time window
%         figure; 
%         set(gcf,'Position',[76 46 411 950]);
%         subplot(311) % LR plots
%         bars = [Delay(dvec(1)).lr_ltw_mn(:,b),Delay(dvec(2)).lr_ltw_mn(:,b),Delay(dvec(3)).lr_ltw_mn(:,b)...
%             Delay(dvec(4)).lr_ltw_mn(:,b),Delay(dvec(5)).lr_ltw_mn(:,b),Delay(dvec(6)).lr_ltw_mn(:,b)];
%         errs = [Delay(dvec(1)).lr_ltw_sem(:,b),Delay(dvec(2)).lr_ltw_sem(:,b),Delay(dvec(3)).lr_ltw_sem(:,b)...
%             Delay(dvec(4)).lr_ltw_sem(:,b),Delay(dvec(5)).lr_ltw_sem(:,b),Delay(dvec(6)).lr_ltw_sem(:,b)];
%        
% %         bars = [Delay(dvec(1)).lr_ltw_mn(evec,b),Delay(dvec(2)).lr_ltw_mn(evec,b),Delay(dvec(3)).lr_ltw_mn(evec,b)...
% %             Delay(dvec(4)).lr_ltw_mn(evec,b),Delay(dvec(5)).lr_ltw_mn(evec,b),Delay(dvec(6)).lr_ltw_mn(evec,b)];
% %         errs = [Delay(dvec(1)).lr_ltw_sem(evec,b),Delay(dvec(2)).lr_ltw_sem(evec,b),Delay(dvec(3)).lr_ltw_sem(evec,b)...
% %             Delay(dvec(4)).lr_ltw_sem(evec,b),Delay(dvec(5)).lr_ltw_sem(evec,b),Delay(dvec(6)).lr_ltw_sem(evec,b)];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         
%         legend('.5','1','2','5','10','20','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('ETW Low Reward Power');
%         hold off; box off
%         %%
%        figure;
%        subplot(211)
%        %HR plot
%         bars = [Delay(dvec(1)).hr_ltw_mn(:,b),Delay(dvec(2)).hr_ltw_mn(:,b),Delay(dvec(3)).hr_ltw_mn(:,b)...
%             Delay(dvec(4)).hr_ltw_mn(:,b),Delay(dvec(5)).hr_ltw_mn(:,b),Delay(dvec(6)).hr_ltw_mn(:,b)];
%         errs = [Delay(dvec(1)).hr_ltw_sem(:,b),Delay(dvec(2)).hr_ltw_sem(:,b),Delay(dvec(3)).hr_ltw_sem(:,b)...
%             Delay(dvec(4)).hr_ltw_sem(:,b),Delay(dvec(5)).hr_ltw_sem(:,b),Delay(dvec(6)).hr_ltw_sem(:,b)];
%      
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
% %        xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('High Reward Power');
%         hold off; box off
%         
%         subplot(212) % HR - LR Difference
%         bars = [Delay(dvec(1)).df_ltw_mn(:,b),Delay(dvec(2)).df_ltw_mn(:,b),Delay(dvec(3)).df_ltw_mn(:,b)...
%             Delay(dvec(4)).df_ltw_mn(:,b),Delay(dvec(5)).df_ltw_mn(:,b),Delay(dvec(6)).df_ltw_mn(:,b)];
%         errs = [Delay(dvec(1)).df_ltw_sem(:,b),Delay(dvec(2)).df_ltw_sem(:,b),Delay(dvec(3)).df_ltw_sem(:,b)...
%             Delay(dvec(4)).df_ltw_sem(:,b),Delay(dvec(5)).df_ltw_sem(:,b),Delay(dvec(6)).df_ltw_sem(:,b)];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('.5','1','2','5','10','20','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('HR-LR Difference Power');
%         hold off; box off
%         flabel=strcat(freqs(b),{' '},...
%             'ETW Bar Plots');
%       %  suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
%         end
%         
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
%         end
        

%     
%     if(plotFigs1 && b==4) % only plot for beta
%         
%         elec = find(evec == target_e); % plot for OFC (can choose different electrode)
%         
%         l_hr = mn_hr_ld(elec,ttvec);
%         l_lr = mn_lr_ld(elec,ttvec);
%         l_df = mn_rewd_ld(elec,ttvec);
%         
%         m_hr = mn_hr_md(elec,ttvec);
%         m_lr = mn_lr_md(elec,ttvec);
%         m_df = mn_rewd_md(elec,ttvec);
%         
%         h_hr = mn_hr_hd(elec,ttvec);
%         h_lr = mn_lr_hd(elec,ttvec);
%         h_df = mn_rewd_hd(elec,ttvec);
%         
%         a_hr = mn_hr_all(elec,ttvec);
%         a_lr = mn_lr_all(elec,ttvec);
%         a_df = mn_rewd_all(elec,ttvec);
%         
%         l_hr_sem = sem_hr_ld(elec,ttvec);
%         l_lr_sem = sem_lr_ld(elec,ttvec);
%         l_df_sem = sem_rewd_ld(elec,ttvec);
%         
%         m_hr_sem = sem_hr_md(elec,ttvec);
%         m_lr_sem = sem_lr_md(elec,ttvec);
%         m_df_sem = sem_rewd_md(elec,ttvec);
%         
%         h_hr_sem = sem_hr_hd(elec,ttvec);
%         h_lr_sem = sem_lr_hd(elec,ttvec);
%         h_df_sem = sem_rewd_hd(elec,ttvec);
%         
%         a_hr_sem = sem_hr_all(elec,ttvec);
%         a_lr_sem = sem_lr_all(elec,ttvec);
%         a_df_sem = sem_rewd_all(elec,ttvec);
%          
%         figure;
%         set(gcf,'Position',[20 13 1610 961]);
%         twlab1 = -500:(3001/length(ttvec)):2500;
%         colvec = ['b','r','g','m'];
%         
%         % shaded err bar plots across 3 delay conditions
%         subplot(341)
%         hold on;
%         shadedErrorBar(twlab1,l_hr(:),l_hr_sem(:),...
%             'lineprops',colvec(1),'patchSaturation',0.33);
%         shadedErrorBar(twlab1,l_lr(:),l_lr_sem(:),...
%             'lineprops',colvec(2),'patchSaturation',0.33); 
%         yl = ylim;
%         plot([0,0,yl(1),yl(2)],'--k')
%         xlabel('Time(ms)');
%         ylabel('Zscored Power'); 
%         title('Low Delay')
%         %legend('HR','LR','Location','best'); 
%         subplot(345)
%         hold on;
%         shadedErrorBar(twlab1,m_hr(:),m_hr_sem(:),...
%             'lineprops',colvec(1),'patchSaturation',0.33);
%         shadedErrorBar(twlab1,m_lr(:),m_lr_sem(:),...
%             'lineprops',colvec(2),'patchSaturation',0.33); 
%         yl = ylim;
%         plot([0,0,yl(1),yl(2)],'--k')
%         xlabel('Time(ms)');
%         ylabel('Zscored Power'); 
%         title('Medium Delay')
%         %legend('HR','LR','Location','best');
%         subplot(349)
%         hold on;
%         shadedErrorBar(twlab1,h_hr(:),h_hr_sem(:),...
%             'lineprops',colvec(1),'patchSaturation',0.33);
%         shadedErrorBar(twlab1,h_lr(:),h_lr_sem(:),...
%             'lineprops',colvec(2),'patchSaturation',0.33);  
%         yl = ylim;
%         plot([0,0,yl(1),yl(2)],'--k')
%         xlabel('Time(ms)');
%         ylabel('Zscored Power'); 
%         title('High Delay')
%         legend('HR','LR','Location','best');
%         
%         % traces plotted for each delay condition(l,m,h) LR, HR, Difference
%         % plot
%         subplot(342) % LR plot
%         plot(tt,l_lr,'Color',[0 .1 .1]);  hold on
%         plot(tt,m_lr,'Color',[0 .4 .4]);  hold on
%         plot(tt,h_lr,'Color',[0 .8 .8]);  hold on;
%         title('Low Reward Traces across Delays');
%         xlabel('Time(ms)');
%         ylabel('Z-score');
%         axis([-500 2500 -2 2]);
%         yl = ylim;
%         plot([0,0,yl(1),yl(2)],'--k')
%         plot([tt(1);tt(end)],[0;0],'k--');
%         box off;                
%         %legend('LD','MD','HD');
%         
%         subplot(343) % HR plot
%         plot(tt,l_hr,'Color',[0 .1 .1]);  hold on
%         plot(tt,m_hr,'Color',[0 .4 .4]);  hold on
%         plot(tt,h_hr,'Color',[0 .8 .8]);  hold on;
%         title('High Reward Traces across Delays');
%         xlabel('Time(ms)');
%         ylabel('Z-score');
%         axis([-500 2500 -2 2]);
%         yl = ylim;
%         plot([0,0,yl(1),yl(2)],'--k')
%         plot([tt(1);tt(end)],[0;0],'k--');
%         box off;                
%         %legend('LD','MD','HD');
%         
%         subplot(344) % HR-LR difference plot
%         plot(tt,l_df,'Color',[0 .1 .1]);  hold on
%         plot(tt,m_df,'Color',[0 .4 .4]);  hold on
%         plot(tt,h_df,'Color',[0 .8 .8]);  hold on;
%         title('High - Low Rew Diff Across Delays');
%         xlabel('Time(ms)');
%         ylabel('Z-score');
%         axis([-500 2500 -2 2]);
%         yl = ylim;
%         plot([0,0,yl(1),yl(2)],'--k')
%         plot([tt(1);tt(end)],[0;0],'k--');
%         box off;                
%         legend('LD','MD','HD','Location','best');
%         
%         %%%%%%
%         
%         % bar plots for individual electrode ETW and LTW
%         % ETW plots - LR 
%         delayList= ["Low","Med","High"];
%         subplot(346)
%         bars =  [mn_l_lr_etw(elec),mn_m_lr_etw(elec),mn_h_lr_etw(elec)]';
%         errs =  [mn_l_lr_etw_sem(elec),mn_m_lr_etw_sem(elec),mn_h_lr_etw_sem(elec)];
%         b1 = bar(bars,'grouped'); hold on;
%         errorbar(bars, errs, 'r' , 'linestyle', 'none');
% %         [ngroups,nbars] = size(bars);
% %         groupwidth = min(0.8,nbars/(nbars+1.5));
% %         for z=1:nbars
% %             % Calculate center of each bar
% %             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
% %             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
% %         end
% %         legend('LD','MD','HD','Location','best');
% %         xticks(1:numel(elecs));
%          xticklabels(delayList);
% %         xtickangle(45);
%         xlabel('Delay Lengths');
%         ylabel('Zscored Mean Power');
%         title('ETW Low Reward Power');
%         hold off; box off
%         
%         subplot(347)
%         bars =  [mn_l_hr_etw(elec),mn_m_hr_etw(elec),mn_h_hr_etw(elec)];
%         errs =  [mn_l_hr_etw_sem(elec),mn_m_hr_etw_sem(elec),mn_h_hr_etw_sem(elec)];
%         b1 = bar(bars,'grouped'); hold on;
%         errorbar(bars, errs, 'r' , 'linestyle', 'none');
% %         [ngroups,nbars] = size(bars);
% %         groupwidth = min(0.8,nbars/(nbars+1.5));
% %         for z=1:nbars
% %             % Calculate center of each bar
% %             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
% %             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
% %         end
% %         legend('LD','MD','HD','Location','best');
% %         xticks(1:numel(elecs));
%         xticklabels(delayList);
% %         xtickangle(45);
%         xlabel('Delay Lengths');
%         ylabel('Zscored Mean Power');
%         title('ETW High Reward Power');
%         hold off; box off
%         
%         subplot(348)
%         bars =  [mn_l_df_etw(elec),mn_m_df_etw(elec),mn_h_df_etw(elec)];
%         errs =  [mn_l_df_etw_sem(elec),mn_m_df_etw_sem(elec),mn_h_df_etw_sem(elec)];
%         b1 = bar(bars,'grouped'); hold on;
%         errorbar(bars, errs, 'r' , 'linestyle', 'none');
% %         [ngroups,nbars] = size(bars);
% %         groupwidth = min(0.8,nbars/(nbars+1.5));
% %         for z=1:nbars
% %             % Calculate center of each bar
% %             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
% %             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
% %         end
% %         legend('LD','MD','HD','Location','best');
% %         xticks(1:numel(elecs));
%         xticklabels(delayList);
% %         xtickangle(45);
%         xlabel('Delay Lengths');
%         ylabel('Zscored Mean Power');
%         title('ETW HR-LR Difference Power');
%         hold off; box off
%         
%         % LTW plots - LR 
%         subplot(3,4,10)
%         bars =  [mn_l_lr_ltw(elec),mn_m_lr_ltw(elec),mn_h_lr_ltw(elec)];
%         errs =  [mn_l_lr_ltw_sem(elec),mn_m_lr_ltw_sem(elec),mn_h_lr_ltw_sem(elec)];
%         b1 = bar(bars,'grouped'); hold on;
%         errorbar(bars, errs, 'r' , 'linestyle', 'none');
% %         [ngroups,nbars] = size(bars);
% %         groupwidth = min(0.8,nbars/(nbars+1.5));
% %         for z=1:nbars
% %             % Calculate center of each bar
% %             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
% %             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
% %         end
% %         legend('LD','MD','HD','Location','best');
% %         xticks(1:numel(elecs));
%         xticklabels(delayList);
% %         xtickangle(45);
%         xlabel('Delay Lengths');
%         ylabel('Zscored Mean Power');
%         title('LTW Low Reward Power');
%         hold off; box off
%         
%         subplot(3,4,11)
%         bars =  [mn_l_hr_ltw(elec),mn_m_hr_ltw(elec),mn_h_hr_ltw(elec)];
%         errs =  [mn_l_hr_ltw_sem(elec),mn_m_hr_ltw_sem(elec),mn_h_hr_ltw_sem(elec)];
%         b1 = bar(bars,'grouped'); hold on;
%         errorbar(bars, errs, 'r' , 'linestyle', 'none');
% %         [ngroups,nbars] = size(bars);
% %         groupwidth = min(0.8,nbars/(nbars+1.5));
% %         for z=1:nbars
% %             % Calculate center of each bar
% %             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
% %             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
% %         end
% %         legend('LD','MD','HD','Location','best');
% %         xticks(1:numel(elecs));
%         xticklabels(delayList);
% %         xtickangle(45);
%         xlabel('Delay Lengths');
%         ylabel('Zscored Mean Power');
%         title('LTW High Reward Power');
%         hold off; box off
%         
%         subplot(3,4,12)
%         bars =  [mn_l_df_ltw(elec),mn_m_df_ltw(elec),mn_h_df_ltw(elec)];
%         errs =  [mn_l_df_ltw_sem(elec),mn_m_df_ltw_sem(elec),mn_h_df_ltw_sem(elec)];
%         b1 = bar(bars,'grouped'); hold on;
%         errorbar(bars, errs, 'r' , 'linestyle', 'none');
% %         [ngroups,nbars] = size(bars);
% %         groupwidth = min(0.8,nbars/(nbars+1.5));
% %         for z=1:nbars
% %             % Calculate center of each bar
% %             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
% %             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
% %         end
% %         legend('LD','MD','HD','Location','best');
% %         xticks(1:numel(elecs));
%         xticklabels(delayList);
% %         xtickangle(45);
%         xlabel('Delay Lengths');
%         ylabel('Zscored Mean Power');
%         title('LTW HR-LR Difference Power');
%         hold off; box off
%         flabel=strcat(freqs(b),{' '},'Electrode',num2str(target_e),{' '},'Reward and Delay Comparisons');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
%         end
%         
%          % Bar plots across all rew related regions
%         figure; 
%         set(gcf,'Position',[76 46 411 950]);
%         subplot(311) % LR plots
%         bars =  [mn_l_lr_etw,mn_m_lr_etw,mn_h_lr_etw];
%         errs =  [mn_l_lr_etw_sem,mn_m_lr_etw_sem,mn_h_lr_etw_sem];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('LD','MD','HD','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('ETW Low Reward Power');
%         hold off; box off
%         
%         subplot(312) % HR bar plot 
%         bars =  [mn_l_hr_etw,mn_m_hr_etw,mn_h_hr_etw];
%         errs =  [mn_l_hr_etw_sem,mn_m_hr_etw_sem,mn_h_hr_etw_sem];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('LD','MD','HD','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('ETW High Reward Power');
%         hold off; box off
%         
%         subplot(313) % HR - LR Difference
%         bars =  [mn_l_df_etw,mn_m_df_etw,mn_h_df_etw];
%         errs =  [mn_l_df_etw_sem,mn_m_df_etw_sem,mn_h_df_etw_sem];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('LD','MD','HD','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('ETW HR-LR Difference Power');
%         hold off; box off
%         flabel=strcat(freqs(b),{' '},...
%             'ETW Bar Plots');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
%         end
%         
%         % Bar plots across all rew related regions LTW
%         figure; 
%         set(gcf,'Position',[76 46 411 950]);
%         subplot(311) % LR plots
%         bars =  [mn_l_lr_ltw,mn_m_lr_ltw,mn_h_lr_ltw];
%         errs =  [mn_l_lr_ltw_sem,mn_m_lr_ltw_sem,mn_h_lr_ltw_sem];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('LD','MD','HD','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('LTW Low Reward Power');
%         hold off; box off
%         
%         subplot(312) % HR bar plot 
%         bars =  [mn_l_hr_ltw,mn_m_hr_ltw,mn_h_hr_ltw];
%         errs =  [mn_l_hr_ltw_sem,mn_m_hr_ltw_sem,mn_h_hr_ltw_sem];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('LD','MD','HD','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('LTW High Reward Power');
%         hold off; box off
%         
%         subplot(313) % HR - LR Difference
%         bars =  [mn_l_df_ltw,mn_m_df_ltw,mn_h_df_ltw];
%         errs =  [mn_l_df_ltw_sem,mn_m_df_ltw_sem,mn_h_df_ltw_sem];
%         b1 = bar(bars,'grouped'); hold on;
% %         errorbar(bars, errs, 'r' , 'linestyle', 'none');
%         [ngroups,nbars] = size(bars);
%         groupwidth = min(0.8,nbars/(nbars+1.5));
%         for z=1:nbars
%             % Calculate center of each bar
%             x = (1:ngroups) - groupwidth/2 + (2*z-1) * groupwidth / (2*nbars);
%             errorbar(x, bars(:,z), errs(:,z), 'r', 'linestyle', 'none');
%         end
%         legend('LD','MD','HD','Location','best');
%         xticks(1:numel(elecs));
%         xticklabels(elecs);
%         xtickangle(45);
%         xlabel('Electrodes');
%         ylabel('Zscored Mean Power');
%         title('LTW HR-LR Difference Power');
%         hold off; box off
%         flabel=strcat(freqs(b),{' '},...
%             'LTW Bar Plots');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir2,char(flabel)),'fig');
%         end
%         
%     end
%     
%     
%     % Plot imagesc plots across all electrodes 
%     if(plotFigs)
%         l_hr = mn_hr_ld(:,ttvec);
%         l_lr = mn_lr_ld(:,ttvec);
%         l_df = mn_rewd_ld(:,ttvec);
%         
%         m_hr = mn_hr_md(:,ttvec);
%         m_lr = mn_lr_md(:,ttvec);
%         m_df = mn_rewd_md(:,ttvec);
%         
%         h_hr = mn_hr_hd(:,ttvec);
%         h_lr = mn_lr_hd(:,ttvec);
%         h_df = mn_rewd_hd(:,ttvec);
%         
%         a_hr = mn_hr_all(:,ttvec);
%         a_lr = mn_lr_all(:,ttvec);
%         a_df = mn_rewd_all(:,ttvec);
% 
%         % Low Delay grouping imagesc plot across REW elecs
%         figure;
%         subplot(311)
%         imagesc(l_hr);colorbar; 
%         title('High Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(312)
%         imagesc(l_lr);colorbar; 
%         title('Low Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(313)
%         imagesc(l_df);colorbar; 
%         title('HR-LR Difference');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         flabel=strcat('Low Delay',freqs(b),{' '},'Mean Power Spectrogram');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir,char(flabel)),'fig');
%         end
%         close;
%         
%         % Medium Delay grouping imagesc plot across REW elecs
%         figure;
%         subplot(311)
%         imagesc(m_hr);colorbar; 
%         title('High Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(312)
%         imagesc(m_lr);colorbar; 
%         title('Low Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(313)
%         imagesc(m_df);colorbar; 
%         title('HR-LR Difference');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         flabel=strcat('Medium Delay',freqs(b),{' '},'Mean Power Spectrogram');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir,char(flabel)),'fig');
%         end
%         close;
%         
%         % High Delay grouping imagesc plot across REW elecs
%         figure;
%         subplot(311)
%         imagesc(h_hr);colorbar; 
%         title('High Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(312)
%         imagesc(h_lr);colorbar; 
%         title('Low Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(313)
%         imagesc(h_df);colorbar; 
%         title('HR-LR Difference');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         flabel=strcat('High Delay',freqs(b),{' '},'Mean Power Spectrogram');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir,char(flabel)),'fig');
%         end
%         close;
%         
%         % All Delays grouping imagesc plot across REW elecs
%         figure;
%         subplot(311)
%         imagesc(h_hr);colorbar; 
%         title('High Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(312)
%         imagesc(h_lr);colorbar; 
%         title('Low Reward');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         subplot(313)
%         imagesc(h_df);colorbar; 
%         title('HR-LR Difference');
%         xlabel('Time');
%         ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%         xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%         xtlbl = linspace(tw1, tw2, numel(xt2));  
%         set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%         set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%         flabel=strcat('All Delays',freqs(b),{' '},'Mean Power Spectrogram');
%         suptitle(flabel);
%         if(saveFigs)
%             saveas(gcf,fullfile(pdir,char(flabel)),'fig');
%         end
%         close;
%         clear l_lr l_hr l_df m_lr m_hr m_df h_lr h_hr h_df a_lr a_hr a_df
%     end
% end
% toc 
% disp('Done grouped imagesc plotting');
% 
% %% For all Delays- Spectogram plots (elecs vs time) for each frequency 
% pdir = 'E:\DD_PhysProcessed\rewL plots\spectAllEs\rewregions_zscored\betaplts\';
% if ~exist(pdir,'dir')
%     mkdir(pdir)
% end
% 
% load('E:\DD_PhysProcessed\times_freqs.mat');
% % Desired Time Window for plotting 
% 
% tw1 = -500;
% tw2 = 2500;
% 
% twa1 = tw1-30;
% twa2 = tw2+30;
% tt = t(t>twa1 & t<twa2);
% ttvec = find(t==min(tt)):find(t==max(tt));
% 
% % Load time and frequency vectors for plotting
% tic
% plotFigs =1;
% evec = [1,5,6,7,25,26,27,28,29];
% freqs = ["Delta","Theta","Alpha","Beta","LGamma","HGamma"];
% elecs = ["M2","A32D","A32V","vOFC","LFC","ALM","LOFC","AIns","VMS","NAcC","NAcS","BLA"];
% twlab1 = tw1:((tw2-tw1+1)/length(ttvec)):tw2;
% evec = [28,27,26,25,31,32,29,30,6,5,1,20]; % 10 target regions in same order as Go/Nogo figure 
% % iterate and plot for signficant electrodes
% 
% for i = 1:length(Delay)-1  % Iterate through all delay lengths 
%     for b = 4 % Fbands 
%         ecnt=0;
%         for e = evec % Electrodes
%             ecnt = ecnt+1;
%             numAni = size(Delay(i).rewL_hr,1);
%             tmp_h = cell2mat(Delay(i).rewL_hr(:,e,b));
%             tmp_l = cell2mat(Delay(i).rewL_lr(:,e,b));
%             plhr(b,ecnt,:) = mean(tmp_h,'omitnan');
%             sehr(b,ecnt,:) = std(tmp_h,'omitnan')/sqrt(numAni);
%             pllr(b,ecnt,:) = mean(tmp_l,'omitnan');
%             selr(b,ecnt,:) = std(tmp_l,'omitnan')/sqrt(numAni);
%         end
%         % Plot imagesc plots across all electrodes
%         if(plotFigs && exist('plhr','var') && exist('pllr','var'))
%             tmp1 = squeeze(plhr(b,:,:));
%             tmp2 = squeeze(pllr(b,:,:));
% 
%             tmp1 = tmp1(:,ttvec);
%             tmp2 = tmp2(:,ttvec);
%             
%             
%             %Correct Plot with bin1 and bin2 for all electrodes
%             figure;
%             subplot(211)
%             imagesc(tmp1);colorbar; 
%             title('High Reward Power');
%             xlabel('Time');
%             ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%             xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%             xtlbl = linspace(tw1, tw2, numel(xt2));  
%             set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%             set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%             subplot(212)
%             imagesc(tmp2);colorbar; 
%             title('Low Reward Power');
%             xlabel('Time');
%             ylabel('Electrodes')
% %             axis([tw1 tw2 -inf inf]);
% %             xtickangle(45);
%             xt2 = [0 100/6 200/6 300/6 400/6 500/6 100];
%             xtlbl = linspace(tw1, tw2, numel(xt2));  
%             set(gca, 'XTick',xt2, 'XTickLabel',xtlbl, 'XTickLabelRotation',45);
%             set(gca, 'YTick',1:numel(evec), 'YTickLabel',elecs, 'YTickLabelRotation',30);
%             flabel=strcat(num2str(Delay(i).length),'s Delay',freqs(b),...
%                 {' '},'Mean Power Spectrogram');
%             suptitle(flabel);
%             if(saveFigs)
%                 saveas(gcf,fullfile(pdir,char(flabel)),'tiff');
%             end
%             close all;
%             clear tmp1 tmp2
%         end
%     end
%     % Reset variables for each delay length
%     clear plhr pllr plhritc pllritc;
%     clear sehr selr sehritc selritc;
% end
% 
% toc
% disp('Done plotting DD rewL electrode spectrograms')