 
function [LFP2,behm,cont]=getmarkers_DD2(sourceFilename,savedir)
% choose initial parameters for reaction time
    try
        [streams,fileh] = load_xdf(sourceFilename);
         cont=0; % if there is some error, cont will equal NaN
    catch
        streams=0;
        behm=0;
        LFP2=0;
        cont=NaN;
        disp('load xdf error');
    end
    if ~isnan(cont)

        RT_highDel = []; RT_lowDel = [];
        while (~isnan(cont))
            if isempty(streams{1,1}.time_series)
                cont=NaN;
                LFP=0;
                LFP2=0;
                behm=0;
                break;
            elseif size(streams,2)<2
                cont=NaN;
                LFP=0;
                LFP2=0;
                Evn=1;
            else
                if length(streams{1,1}.time_series(:,1))<2 % specify which stream is LFP and which one is events
                    LFP=2;
                    Evn=1;
                else
                    LFP=1;
                    Evn=2;
                end

                if  streams{1,LFP}.info.effective_srate>1200   % get rid of unknown sample rate
                    cont=NaN;
                    LFP2=0;
                    behm=0;
                    break;
                end

                if isempty(streams{1,LFP}.time_stamps) || isempty(streams{1,Evn}.time_stamps)
                    cont=NaN;
                    LFP2=0;
                    behm=0;
                    break;
                end
                if size(streams,2)>=2
                    %change sec to ms to match with eeglab              
                    streams{1,Evn}.time_stamps=(streams{1,Evn}.time_stamps-streams{1,LFP}.time_stamps(1))*1000;
                    streams{1,LFP}.time_stamps=(streams{1,LFP}.time_stamps-streams{1,LFP}.time_stamps(1))*1000;

                end
            end

            %% new DD data
            tmp=streams{1,Evn}.time_series;
            tmp2=streams{1,Evn}.time_stamps;

            % Assign LSL event codes

            sind=find(tmp==12);
            sindT=tmp2(sind);

            sind2=find(tmp==14);
            sindT2=tmp2(sind2);

            Rew3=(find(tmp==3));
            rew3T=tmp2(Rew3);


            %%
            trial = zeros(length(sind),1);
            reentry = zeros(length(sind),1);
            dec = zeros(length(sind),1);
            rTime = zeros(length(sind),1);
            skptr = zeros(length(sind),1);
            swap = zeros(length(sind),1);
            swapT=zeros(length(sind),1);
            delayDur2 = zeros(length(sind),1); 
            delayDur4 = zeros(length(sind),1);
            count2 = 0;
            count4 = 0;
            entry =0;
            %%
            for i = 5:length(sind)
                entry = 1;
                %respTime = Rew_time(i)-sindT(i);
                if i<length(sind)
                    tmp3=tmp([sind(i):sind(i+1)-1]);
                    tmp4=tmp2([sind(i):sind(i+1)-1]);
                else
                    tmp3=tmp([sind(i):end]);
                    tmp4=tmp2([sind(i):end]);
                end

                %IR sensors
                NP2=find(tmp3==102);
                NP3=find(tmp3==103);
                NP4=find(tmp3==104);

                % Reward timing
                Rew3=find(tmp3==3);

                %IR sensors timing
                NP3Time = tmp4(NP3);
                NP4Time = tmp4(NP4);
                NP2Time = tmp4(NP2);

                if ~isempty(Rew3)
                    Rew3Time(i) = tmp4(Rew3(1));
                    tmprr=find(NP3>Rew3(1));
                    if ~isempty(tmprr)
                        RewrespT(i)=tmp4(NP3(tmprr(1)));%-tmp4(1);
                        %RewrespT(i) = tmp4(NP3(tmprr(1)))- Rew3Time(i);
                    else
                        RewrespT(i)=NaN;
                    end
                else
                    Rew3Time(i) = NaN;
                    RewrespT(i)=NaN;      
                end


                %after 60s no response-- trial skipped
                skptr(i) = 0;
                check2 = 0;
                if(isempty(NP2) && isempty(NP4))
                    skptr(i) = 1;
                    trial(i) = 1;
                    Rew3Time(i) = NaN;

                elseif(isempty(NP2))
                    trial(i) = 4;
                    count4 = count4 +1;
                elseif(isempty(NP4))
                    trial(i) = 2;
                    count2 = count2 + 1;
                elseif((NP4Time(1) < NP2Time(1)))
                    trial(i) = 4;
                    count4 = count4 +1;
                    check2 = 1; %Both NP2 and NP4 become activated
                elseif(NP2Time(1) < NP4Time(1))
                    trial(i) = 2;
                    count2 = count2 + 1;
                    check2 = 1;
                end

                res_2{i}=NP2Time;
                res_4{i}=NP4Time;

                % NP4 Trial(High Delay)
                if(trial(i)==4)
                    rTime(i) = NP4Time(1)-tmp4(1);
                    delayDur4(i) = Rew3Time(i)-NP4Time(1); 
                    if(check2)
                        swapT(i) = NP2Time(1)-tmp4(1);
                        swap(i) = 1;
                    end
                    %NP2 Trial (Low Delay)
                elseif(trial(i)==2) 
                    rTime(i) = NP2Time(1) - tmp4(1);
                    delayDur2(i) = Rew3Time(i)-NP2Time(1); 
                    if(check2)
                        swapT(i) = NP4Time(1)-tmp4(1);
                        swap(i) = 1;
                    end
                end
            end
            %%
            delayDur2(delayDur2==0) = []; 
            delayDur4(delayDur4==0) = [];

            trType = round((nanmean(delayDur4))); % High delay length 
            if entry
                %keyboard
            behm.resp = trial; %Trial = 2,4 or 1 if tr skipped
            behm.startRT = rTime; %reaction time(time for response after tr start)

            behm.onsetT = sindT; %Trial onset determined by LED2+LED4 being on
            behm.respT = sindT+rTime; %Response onset time
            behm.rewT = Rew3Time; %Reward onset time
            behm.rewrespT = RewrespT; %Reward Response onset time (rew onset to first response)
            behm.rewRT = RewrespT-behm.rewT; %reaction time(from rew to NP3 entry).

            behm.skptr = skptr; % 1 if tr skipped, 0 if normal tr
            behm.swap = swap; %Whether reentry was to the other NP
            behm.swapT = swapT; %Time of response if swapped NPs
            behm.res_2=res_2; % all response times on NP2 on every trial
            behm.res_4=res_4; % all response times on NP4 on every trial

            behm.ld=nanmean(delayDur2); % calculate low reward delay (should be either .5 or 5)
            behm.hd=nanmean(delayDur4);% calculate high reward delay (should be either .5 to 20)

            behm.numHighTr = count4;
            behm.numLowTr = count2;
            behm.trType = trType;
            %clear streams
            if ~exist(savedir,'dir')
                mkdir(savedir);
            end
            if LFP~=0
                LFP2.ts=streams{1,LFP}.time_stamps;
                LFP2.ch=streams{1,LFP}.time_series;
                LFP2.sr=streams{1,LFP}.segments.effective_srate;
            end
            sfile= [savedir 'extracted_data.mat'];
            sfig= [savedir 'beh.tif'];
            figure;
            hist(behm.resp);
            d=hist(behm.resp);
            hold on;
            text(4,d(end)+2,num2str(d(end)));
            text(2,d(5)+2,num2str(d(5)));
            saveas(gcf,sfig,'tiff');
            close
            end
            if exist('behm','var')
                if exist('LFP2','var')
                    save(sfile, 'LFP2','behm','-v7.3');
                else
                    LFP2=0;
                    save(sfile,'behm','-v7.3');
                end
                break;
            else
                cont=NaN;
                behm=0;
                LFP2=0;
            end
        end
    end
end
