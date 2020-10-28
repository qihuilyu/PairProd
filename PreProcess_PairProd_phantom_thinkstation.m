clear
close all
clc

patientName = 'newphantom_10k_newraytrace';
projectName = 'PairProd';
patFolder = fullfile('D:\datatest\PairProd\',patientName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
h5file = fullfile(dosecalcFolder,'PairProd_beamletdose.h5');
maskfile = fullfile(dosecalcFolder,'PairProd_masks.h5');
fmapsfile = fullfile(dosecalcFolder,'PairProd_fmaps.h5');
Anni3Dfile = fullfile(dosecalcFolder,'PairProd_NofPositronAnni3D.h5');
DetectedEventsfile = fullfile(dosecalcFolder,'PairProd_DetectedEvents.h5');

%% masks, fmaps, dose matrix, annihilation matrix
[M,M_Anni,dose_data,masks]=BuildDoseMatrix_PairProd(h5file, maskfile, fmapsfile, Anni3Dfile);  

dose = reshape(full(sum(M,2)),size(masks{1}.mask));
Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));

figure;imshow3D(dose)
figure;imshow3D(Anni3D)

pdose = 25;
[StructureInfo, params] = InitIMRTparams_DLMCforRyan(M,dose_data,masks,pdose,[1,2,3]);

projectFolder = fullfile(patFolder,'PairProd');
paramsFolder = fullfile(projectFolder,'params');
mkdir(paramsFolder)

ParamsNum = 0;
save(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
InfoNum = 0;
save(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');

dosematrixFolder = fullfile(projectFolder,'dosematrix');
mkdir(dosematrixFolder)
save(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks','-v7.3');

%% Dicom
DicomPath = fullfile(dosecalcFolder,'ctdata');
baseFileNames = dir([DicomPath '/*.dcm']);
[sortedFile,index] = sort_nat({baseFileNames.name});

img = [];
for ii = 1:numel(sortedFile)
    dcinfo = dicominfo(fullfile(DicomPath,sortedFile{ii}));
    if(~strcmp(dcinfo.Modality,'RTstruct'))
        imgres = dcinfo.SliceThickness;
        img(:,:,ii) = dicomread(fullfile(DicomPath,sortedFile{ii}));
    end
end
figure;imshow3D(img,[0,2000])
save(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');


%% Ring detection
info = h5info(DetectedEventsfile);
detectorIds = double(h5read(DetectedEventsfile, '/detectorIds')) + 1; %
beamNo = double(h5read(DetectedEventsfile, '/beamNo')) + 1; %
beamletNo = double(h5read(DetectedEventsfile, '/beamletNo')) + 1; %
energy = h5read(DetectedEventsfile, '/energy'); %
eventIds = double(h5read(DetectedEventsfile, '/eventIds')) + 1; %
globalTimes = h5read(DetectedEventsfile, '/globalTimes'); %
Mega = [globalTimes,eventIds,energy,beamletNo,beamNo,detectorIds];
save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection_original.mat']),'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes','-v7.3');

%% Clean data
numevent = max(eventIds);
beamSizes = squeeze(sum(sum(params.BeamletLog0,1),2));
cumsumbeamSizes = cumsum([0; beamSizes]);
beamNoshift = cumsumbeamSizes(beamNo);
beamletIDs = double(beamletNo) + beamNoshift;
AlleventID = (beamletIDs-1)*numevent + eventIds;
nb_cryst = max(detectorIds);
AlldetectorID = (AlleventID-1)*nb_cryst + detectorIds;
[sortedAlldetectorID, sortAlldetectorIDInd] = sort(AlldetectorID);
sortInd_sameparticle = find(diff(sortedAlldetectorID)==1);
Ind_coin1 = sortAlldetectorIDInd(sortInd_sameparticle);
Ind_coin2 = sortAlldetectorIDInd(sortInd_sameparticle+1);
mask_sameenergy = (energy(Ind_coin1)-energy(Ind_coin2)==0);
timediff = globalTimes(Ind_coin1)-globalTimes(Ind_coin2);
badIDbuff1 = Ind_coin1(mask_sameenergy & timediff>0);
badIDbuff2 = Ind_coin2(mask_sameenergy & timediff<0);
badID1 = union(badIDbuff1,badIDbuff2);

AlldetectorID2 = (AlleventID-1)*nb_cryst + mod(detectorIds,nb_cryst);
[sortedAlldetectorID2, sortAlldetectorIDInd2] = sort(AlldetectorID2);
sortInd_sameparticle2 = find(diff(sortedAlldetectorID2)==1);
Ind_coin21 = sortAlldetectorIDInd2(sortInd_sameparticle2);
Ind_coin22 = sortAlldetectorIDInd2(sortInd_sameparticle2+1);
mask_sameenergy = (energy(Ind_coin21)-energy(Ind_coin22)==0);
timediff = globalTimes(Ind_coin21)-globalTimes(Ind_coin22);
badIDbuff1 = Ind_coin21(mask_sameenergy & timediff>0);
badIDbuff2 = Ind_coin22(mask_sameenergy & timediff<0);
badID2 = union(badIDbuff1,badIDbuff2);

badID = union(badID1,badID2);

beamletIDs(badID,:) = [];
detectorIds(badID,:) = [];
beamNo(badID,:) = [];
beamletNo(badID,:) = [];
energy(badID,:) = [];
eventIds(badID,:) = [];
globalTimes(badID,:) = [];
Mega(badID,:) = [];
save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes','-v7.3');

%% fluence map segments
numbeamlets = size(M,2);

x_ = ones(numbeamlets,1);
dose = reshape(M*x_,size(masks{1}.mask));
figure;imshow3D(dose,[])

selectbeamNo = 1;
x_onebeam = zeros(numbeamlets,1);
x_onebeam(cumsumbeamSizes(selectbeamNo)+1:cumsumbeamSizes(selectbeamNo+1)) = 1;
dose_onebeam = reshape(M*x_onebeam,size(masks{1}.mask));
figure;imshow3D(dose_onebeam,[])

%% Time correction: Adding time of previous events
doserate = 0.1/60; % (0.1Gy/min)
time = max(dose_onebeam(:))*numevent/doserate;
eventrate = time/numevent*1e+09; % ns

numeventbatch = 1e+05;
numbatches = numevent/numeventbatch;
deltatime_event = normrnd(eventrate,eventrate/5,numeventbatch,numbeamlets);
cumsum_eventtime_batch = cumsum(deltatime_event);

event_perbatchIDs = mod(eventIds-1, numeventbatch)+1;
event_batchID = (eventIds - event_perbatchIDs)/numeventbatch + 1;
event_perbatch_beamletIDs = sub2ind([numeventbatch,numbeamlets], event_perbatchIDs, beamletIDs);

numbeams = max(beamNo(:));
batchtime = max(cumsum_eventtime_batch(event_perbatch_beamletIDs));
beamtime = batchtime*numbatches;
deltatime_beam = beamtime*(1:numbeams)'; % 1 ms

CorrectedTime = globalTimes + cumsum_eventtime_batch(event_perbatch_beamletIDs)...
        + event_batchID*batchtime + beamNo*beamtime;
[sortedtime, sortInd] = sort(CorrectedTime);
save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'CorrectedTime','sortedtime','sortInd','-append');

% %% Identify LOR
% Ind_511 = find(abs(energy-0.511)<0.0001);
% EnergyResolution = 0.1;
% Ind_accept = find(abs(energy-0.511)<0.511*EnergyResolution);
% 
% CoincidenceTime = 2;  % ns  
% % [sortedtime, sortInd] = sort(CorrectedTime);
% sortInd_coin = find(diff(sortedtime)<CoincidenceTime);
% Ind_coin1 = sortInd(sortInd_coin);
% Ind_coin2 = sortInd(sortInd_coin+1);
% 
% [Ind_coin1_511buff, iInd_coin1, iInd_511_1] = intersect(Ind_coin1, Ind_511);
% [Ind_coin2_511buff, iInd_coin2, iInd_511_2] = intersect(Ind_coin2, Ind_511);
% 
% iInd_coin = intersect(iInd_coin1,iInd_coin2);
% Ind_coin1_511 = Ind_coin1(iInd_coin);
% Ind_coin2_511 = Ind_coin2(iInd_coin);
% 
% [Ind_coin1_acceptbuff, iInd_coin1, iInd_accept_1] = intersect(Ind_coin1, Ind_accept);
% [Ind_coin2_acceptbuff, iInd_coin2, iInd_accept_2] = intersect(Ind_coin2, Ind_accept);
% 
% iInd_coin = intersect(iInd_coin1,iInd_coin2);
% Ind_coin1_accept = Ind_coin1(iInd_coin);
% Ind_coin2_accept = Ind_coin2(iInd_coin);
% 
% TruePositive = length(Ind_coin1_511)/length(Ind_coin1_accept);
% 
% %% Image Reconstruction
% R1 = 120;
% distrange = 30;
% 
% detid_pair = [detectorIds(Ind_coin1_accept) detectorIds(Ind_coin2_accept)];
% Sino = rebin_PET(detid_pair, nb_cryst, R1, distrange);
% figure;imshow(Sino,[])
% 
% ig = image_geom('nx', 151, 'ny', 153, 'fov', 45);
% sg = sino_geom('par', 'nb', size(Sino,2), 'na', nb_cryst, ...
%     'dr', 2*distrange/size(Sino,2));
% img_fbp = em_fbp(sg, ig, Sino');
% figure;imshow(img_fbp/max(img_fbp(:)),[])
% 
% Anni2D = Anni3D(:,:,ceil(end/2));
% figure;imshow([img_fbp/max(img_fbp(:)) Anni2D/max(Anni2D(:))],[])
% 
% 
% %% Reconstruction-less image generation
% 
% imgsize = size(img_fbp);
% reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
%     'imgres',imgres,'imgsize',imgsize);
% 
% Ind_coin = [Ind_coin1_accept Ind_coin2_accept];
% img_direct = PETrecon_TOF_direct(CorrectedTime, detectorIds, Ind_coin, reconparams);
% figure; imshow(img_direct,[])
% figure;imshow([img_fbp/max(img_fbp(:)) Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:))],[])
% 
% 
% % deltat = CorrectedTime(Ind_coin1_accept) - CorrectedTime(Ind_coin2_accept);
% % clight = 30; % c = 30cm/ns
% % deltar = deltat*clight/2; 
% % % figure; hist(deltar)
% % % figure; hist(deltat)
% % 
% % theta = detid_pair/nb_cryst*2*pi;
% % x_ = R1*sin(theta);
% % y_ = R1*cos(theta);
% % 
% % p1 = [x_(:,1), y_(:,1)];
% % p2 = [x_(:,2), y_(:,2)];
% % p1_p2 = p1 - p2;
% % n_p1_p2 = p1_p2./sqrt(sum(p1_p2.^2,2));
% % p0 = -n_p1_p2.*deltar + (p1+p2)/2;
% % 
% % test_all = [detid_pair,p1,p2,p0,sqrt(sum(p0.^2,2))];
% % % figure;hist(sqrt(sum(p0.^2,2)),100)
% % 
% % %% Show image
% % p0(sqrt(sum(p0.^2,2))>80) = NaN;
% % p0(find(isnan(p0(:,1))),:)=[];
% % xp0 = p0(:,1);
% % yp0 = p0(:,2);
% % % figure;hist(sqrt(sum(p0.^2,2)),100)
% % % figure;hist3(p0,'Nbins',[100,100])
% % % figure;scatter(xp0, yp0, 1); axis equal
% % % 
% % % xImage = floor(xp0 - min(xp0)) + 1;
% % % yImage = floor(yp0 - min(yp0)) + 1;
% % 
% % res = 0.3;
% % img_norecon = zeros(size(img,1),size(img,2));
% % xImage = ((1:size(img,1))-(1+size(img,1))/2)*res;
% % yImage = ((1:size(img,2))-(1+size(img,2))/2)*res;
% % % Image = zeros(450,450);
% % % xImage = ((1:450)-(1+450)/2)*res;
% % % yImage = ((1:450)-(1+450)/2)*res;
% % 
% % for ii = 1:length(xp0)
% %     ixp0 = xp0(ii);
% %     iyp0 = yp0(ii);
% %     
% %     if(ixp0<min(xImage)||ixp0>max(xImage)||iyp0<min(yImage)||iyp0>max(yImage))
% %          continue
% %     end
% %     
% %     [~,xind] = min(abs(ixp0-xImage));
% %     [~,yind] = min(abs(iyp0-yImage));
% % 
% %     img_norecon(xind,yind) = img_norecon(xind,yind) + 1;
% % end
% % 
% % figure; imshow(img_norecon,[])
% % figure;imshow([img_fbp/max(img_fbp(:)) Anni2D/max(Anni2D(:)) img_norecon/max(img_norecon(:))],[])
% 
% 
% % %%
% % 
% % [sortedenergy, sortenergyInd] = sort(energy);
% % sortInd_coinenergy = find(diff(sortedenergy)==0);
% % Ind_coin1 = sortenergyInd(sortInd_coinenergy);
% % Ind_coin2 = sortenergyInd(sortInd_coinenergy+1);
% % 
% % 
% % test = [AlleventID(Ind_coin1) AlleventID(Ind_coin2)];
% % test2 = test(1:100,:);
% % 
% % 
% % 
% % 
% 
% 
% 
% 
% 
% 
