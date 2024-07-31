function extractLFPbeh_DD2(anf,savedir,par,rerun)

rdir=anf;  % point to directory where raw data is kept.

%Assuming current organization of raw data is based on data acquiring. Will
%resort/save data for each animal based on data acquired.

sdir=savedir;    %point to directory where processed data will be kept.
if ~exist(sdir)
    mkdir(sdir)
end
if rerun.beh==1
       [LFP,behm,cont]=getmarkers_DD2(anf,savedir); % gets behavioral markers for DD paradigm
    if(cont==0 && size(LFP.ch,1)>30)
        if rerun.ERP==1
        epochdata=LFP_ERP_DD2(LFP,behm,savedir,par);
            if rerun.ERSP==1
                clear LFP behm;
                ERSP_DD2(epochdata,savedir);
            end
        end
    end
        
end
    
end
