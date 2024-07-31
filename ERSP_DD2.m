function ERSP_DD2(behLFP,savedir)
plt=1;
delete(gcp('nocreate'))
parpool(16);

Frame=7001;
Epoches=[-3000,4000];
channelN=32;
Base=[-2500,-1500];
Base1=NaN;
MinF=2; MaxF=150;
SampleRate=1000;
FN={'De', 'Th','Al','Be','G1','G2'};

tmp1=behLFP.tst;
tmp1a=(behLFP.tst(:,:,behLFP.beh.resp==2));
tmp1b=(behLFP.tst(:,:,behLFP.beh.resp==4));
bd1=squeeze(isnan(tmp1a(25,3000,:)));
bd2=squeeze(isnan(tmp1b(25,3000,:)));
tmp1a(:,:,bd1)=[];
tmp1b(:,:,bd2)=[];

tmp2=behLFP.resp;
tmp2a=behLFP.resp(:,:,behLFP.beh.resp==2);
tmp2b=(behLFP.resp(:,:,behLFP.beh.resp==4));
bd1=squeeze(isnan(tmp2a(25,3000,:)));
bd2=squeeze(isnan(tmp2b(25,3000,:)));
tmp2a(:,:,bd1)=[];
tmp2b(:,:,bd2)=[];


tmp3=behLFP.rew;
tmp3a=(behLFP.rew(:,:,behLFP.beh.resp==2));
tmp3b=(behLFP.rew(:,:,behLFP.beh.resp==4));
bd1=squeeze(isnan(tmp3a(25,3000,:)));
bd2=squeeze(isnan(tmp3b(25,3000,:)));
tmp3a(:,:,bd1)=[];
tmp3b(:,:,bd2)=[];

tmp4=behLFP.rewresp; 
%tmp4a=(behLFP.rewresp(:,:,behLFP.rewresp==2)); 
%tmp4b=(behLFP.rewresp(:,:,behLFP.rewresp==4));


%clear behLFP;
if exist('tmp1','var')
    if(size(tmp1,3)>10)
        parfor CH=1:32
            if ~isnan((squeeze(tmp1(CH,300,50))))
                disp(['this elect is running ', num2str(CH)]);
                [~,~,~,~,~,~,~,DD_st(CH,:,:,:)]=newtimef(squeeze(tmp1(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base1, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [~,~,~,~,~,~,~,DD_resp(CH,:,:,:)]=newtimef(squeeze(tmp2(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base1, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [~,~,~,t(CH,:),freqs(CH,:),~,~,DD_rew(CH,:,:,:)]=newtimef(squeeze(tmp3(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base1, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [~,~,~,t(CH,:),freqs(CH,:),~,~,DD_rewresp(CH,:,:,:)]=newtimef(squeeze(tmp4(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base1, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                
                [DD_lrstp(CH,:,:),DD_lrstitc(CH,:,:),~,~,~,~,~,~]=newtimef(squeeze(tmp1a(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [DD_lrresp(CH,:,:),DD_lrresitc(CH,:,:),~,~,~,~,~,~]=newtimef(squeeze(tmp2a(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [DD_lrrewp(CH,:,:),DD_lrrewitc(CH,:,:),~,~,~,~,~,~]=newtimef(squeeze(tmp3a(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                
                [DD_hrstp(CH,:,:),DD_hrstitc(CH,:,:),~,~,~,~,~,~]=newtimef(squeeze(tmp1b(CH,:,:)), Frame, Epoches,  SampleRate, [2 0.7],'baseline', Base, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [DD_hrresp(CH,:,:),DD_hrresitc(CH,:,:),~,~,~,~,~,~]=newtimef(squeeze(tmp2b(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
                [DD_hrrewp(CH,:,:),DD_hrrewitc(CH,:,:),~,~,~,~,~,~]=newtimef(squeeze(tmp3b(CH,:,:)), Frame, Epoches, SampleRate, [2 0.7],'baseline', Base, 'freqs', [2 150],'plotersp','off','plotitc','off','verbose','off');
            end
        end
        freqs=freqs(1,:);
        t=t(1,:);

        save([savedir 'DD_st.mat'], 'DD_st','DD_hrstp','DD_hrstitc','DD_lrstp','DD_lrstitc','t','freqs','-v7.3');
        save([savedir 'DD_resp.mat'], 'DD_resp','DD_hrresp','DD_hrresitc','DD_lrresp','DD_lrresitc','t','freqs','-v7.3');
        save([savedir 'DD_rew.mat'], 'DD_rew','DD_hrrewp','DD_hrrewitc','DD_lrrewp','DD_lrrewitc','t','freqs','-v7.3' );
        save([savedir 'DD_rewresp.mat'], 'DD_rewresp','t','freqs','-v7.3' );

        if plt==1
                F{1}=find(freqs>1 & freqs<4);
                F{2}=find(freqs>3 & freqs<7);
                F{3}=find(freqs>8 & freqs<12);
                F{4}=find(freqs>15 & freqs<30);
                F{5}=find(freqs>40 & freqs<65);
                F{6}=find(freqs>80 & freqs<150);
                
                for i = 1:length(F)
                    tmp1e=squeeze(mean((DD_lrstp(:,F{i},:)),2));
                    tmp3e=squeeze(mean((DD_hrstp(:,F{i},:)),2));
                    tmp1i=squeeze(mean(abs(DD_lrstitc(:,F{i},:)),2));
                    tmp3i=squeeze(mean(abs(DD_hrstitc(:,F{i},:)),2));

                    try
                    figure
                    subplot(221)
                    imagesc(t,[],tmp1e); colorbar; caxis([-5 5]); title([FN{i} ' LR Trial Lock']);
                    subplot(223)
                    imagesc(t,[],tmp3e); colorbar;  caxis([-5 5]); title([FN{i} 'HR Trial Lock']);
                    subplot(222)
                    imagesc(t,[],tmp1i); colorbar; caxis([-.5 .5]);
                    subplot(224)
                    imagesc(t,[],tmp3i); colorbar; caxis([-.5 .5]);
                    saveas(gcf,[savedir 'plots/' FN{i} ' DD_TL.tif'],'tif')
                    
                    tmp1e=squeeze(mean((DD_lrresp(:,F{i},:)),2));
                    tmp3e=squeeze(mean((DD_hrresp(:,F{i},:)),2));
                    tmp1i=squeeze(mean(abs(DD_lrresitc(:,F{i},:)),2));
                    tmp3i=squeeze(mean(abs(DD_hrresitc(:,F{i},:)),2));
                    
                    figure
                    subplot(221)
                    imagesc(t,[],tmp1e); colorbar;  caxis([-5 5]); title([FN{i} ' LR Resp Lock']);
                    subplot(223)
                    imagesc(t,[],tmp3e); colorbar; caxis([-5 5]);  title([FN{i} 'HR Resp Lock']);
                    subplot(222)
                    imagesc(t,[],tmp1i); colorbar; caxis([-.5 .5]);
                    subplot(224)
                    imagesc(t,[],tmp3i); colorbar; caxis([-.5 .5]);
                    saveas(gcf,[savedir 'plots/' FN{i} ' DD_RespL.tif'],'tif')
                    
                    tmp1e=squeeze(mean((DD_lrrewp(:,F{i},:)),2));
                    tmp3e=squeeze(mean((DD_hrrewp(:,F{i},:)),2));
                    tmp1i=squeeze(mean(abs(DD_lrrewitc(:,F{i},:)),2));
                    tmp3i=squeeze(mean(abs(DD_hrrewitc(:,F{i},:)),2));
                    
                    figure
                    subplot(221)
                    imagesc(t,[],tmp1e); colorbar; caxis([-5 5]); title([FN{i} ' LR Rew Lock']);
                    subplot(223)
                    imagesc(t,[],tmp3e); colorbar; caxis([-5 5]); title([FN{i} 'HR Rew Lock']);
                    subplot(222)
                    imagesc(t,[],tmp1i); colorbar; caxis([-.5 .5]);
                    subplot(224)
                    imagesc(t,[],tmp3i); colorbar; caxis([-.5 .5]);
                    saveas(gcf,[savedir 'plots/' FN{i} ' DD_RewL.tif'],'tif')
                    end
                end               
            end
        
        clear D*
        close all;
    end
end

end
