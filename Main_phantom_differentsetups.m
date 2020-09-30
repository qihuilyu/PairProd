% clear
% close all
% clc
% 
% patientName = 'phantom';
% projectName = 'PairProd';
% patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
% projectFolder = fullfile(patFolder,projectName);
% dosecalcFolder = fullfile(patFolder,'dosecalc');
% dosematrixFolder = fullfile(projectFolder,'dosematrix');
% resultFolder = fullfile(projectFolder,'result');
% mkdir(resultFolder)
% 
% load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']));
% CorrectedTime0 = CorrectedTime;
% detectorIds0 = detectorIds;
% energy0 = energy;
% numevents0 = max(eventIds);
% load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
% load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');
% 
% paramsFolder = fullfile(projectFolder,'params');
% ParamsNum = 0;
% load(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
% InfoNum = 0;
% load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');
% 
% mask = StructureInfo(2).Mask(:,:,56);
% mumap = double(mask);
% mumap(mask==1) = 0.0096;
% mumap(mask==0) = 0.0001;
% 
% dose3D = reshape(full(sum(M,2)),size(masks{1}.mask));
% dose2D = dose3D(:,:,ceil(end/2));
% 
% Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
% Anni2D = Anni3D(:,:,ceil(end/2));

%% Identify LOR
EnergyResolution = 0.1;
CoincidenceTime = 2;  % ns

% flag = 0;
profile on
for DS = [1,0.1]
    numevents = DS*numevents0;
    indDS = find(eventIds<=numevents);
    
    CorrectedTimebuff = CorrectedTime0(indDS);
    detectorIds = detectorIds0(indDS);
    energy = energy0(indDS);
    
    for TR = [0, 0.02, 0.1]
        CorrectedTime = CorrectedTimebuff + TR*randn(size(CorrectedTimebuff));
        [sortedtime, sortInd] = sort(CorrectedTime);
        
        Ind_coin_511 = IdentifyLOR_511(energy, sortedtime, sortInd, CoincidenceTime);
        Ind_coin_accept = IdentifyLOR(energy, sortedtime, sortInd, CoincidenceTime, EnergyResolution);
        TruePositive = length(Ind_coin_511)/length(Ind_coin_accept);
        
        detid_pair = detectorIds(Ind_coin_accept);
        R1 = 1200;
        distrange = 300;
        imgsize = size(img);
        nb_cryst = max(detectorIds);
        [sino, dr, sinobuff] = rebin_PET(detid_pair, nb_cryst, R1, distrange);
        if(flag==0)
            ig = image_geom('nx', imgsize(1), 'ny', imgsize(2), 'fov', imgsize(1)*imgres);
            sg = sino_geom('par', 'nb', size(sino,1), 'na', nb_cryst, 'dr', dr);
            figure; sg.plot([ig])
            G = Gtomo2_strip(sg, ig);
            ci = GetACfactor_sino(G, mumap);
            
            flag=1;
        end
        figure;imshow(sino,[])
        saveas(gcf,fullfile(resultFolder,['sino_TR_' num2str(TR) '_nevents_' num2str(numevents) '.png']))
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        imshow(dose2D*numevents*10^3,[]); colorbar; title(['Dose (mGy): maximum at ' num2str(max(dose2D(:))*numevents*10^3) ' mGy'])
        set(gcf,'outerposition',[0 0 1 1]);
        saveas(gcf,fullfile(resultFolder,['dose2D_TR_' num2str(TR) '_nevents_' num2str(numevents) '.png']))
        
        img_fbp = em_fbp_QL(sg, ig, sino./ci);
        cilist = GetACfactor_list(ci, detid_pair, nb_cryst, R1, distrange);
        reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
            'imgres',imgres,'imgsize',size(img_fbp));
        img_direct = recon_TOF_direct(CorrectedTime, detectorIds, Ind_coin_accept, cilist, reconparams);
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        imshow([Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:)) img_fbp/max(img_fbp(:)) ],[])
        set(gcf,'outerposition',[0 0 1 1]);
        saveas(gcf,fullfile(resultFolder,['reconimg_TR_' num2str(TR) '_nevents_' num2str(numevents) '.png']))
        
    end
    
end


profile report
