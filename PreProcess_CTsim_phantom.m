clear
close all
clc
% 
% t = tic;
% while toc(t)<7200*3
%     pause(2)
% end
% !('/media/raid1/qlyu/PairProd/datatest/collect_doescalc_pairprod.sh')

patientName = 'CTphantom_360beam_10m_CTsim';
projectName = 'CTsim';
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
h5file = fullfile(dosecalcFolder,[projectName '_beamletdose.h5']);
maskfile = fullfile(dosecalcFolder,[projectName '_masks.h5']);
fmapsfile = fullfile(dosecalcFolder,[projectName '_fmaps.h5']);
DetectedEventsCTfile = fullfile(dosecalcFolder,[projectName '_DetectedEventsCT.h5']);

%% masks, fmaps, dose matrix, annihilation matrix
[M,dose_data,masks]=BuildDoseMatrix_CTsim(h5file, maskfile, fmapsfile);

dose = reshape(full(sum(M,2)),size(masks{1}.mask));
figure;imshow3D(dose)

pdose = 25;
[StructureInfo, params] = InitIMRTparams_DLMCforRyan(M,dose_data,masks,pdose,[1,2,0]);

projectFolder = fullfile(patFolder,'PairProd');
paramsFolder = fullfile(projectFolder,'params');
mkdir(paramsFolder)

ParamsNum = 0;
save(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
InfoNum = 0;
save(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');

dosematrixFolder = fullfile(projectFolder,'dosematrix');
mkdir(dosematrixFolder)
save(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','dose_data','masks','-v7.3');

%% Dicom
DicomPath = fullfile(dosecalcFolder,'ctdata');
baseFileNames = dir([DicomPath '/*.dcm']);
[sortedFile,index] = sort_nat({baseFileNames.name});

img = [];
for ii = 1:numel(sortedFile)
    dcinfo = dicominfo(fullfile(DicomPath,sortedFile{ii}));
    if(~strcmpi(dcinfo.Modality,'RTstruct'))
        imgres = dcinfo.SliceThickness;
        img(:,:,ii) = dicomread(fullfile(DicomPath,sortedFile{ii}));
    end
end
figure;imshow3D(img,[0,2000])
save(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');


%% Ring detection
info = h5info(DetectedEventsCTfile);
detectorIds1 = double(h5read(DetectedEventsCTfile, '/detectorIds1')) + 1; %
detectorIds2 = double(h5read(DetectedEventsCTfile, '/detectorIds2')) + 1; %
beamNo = double(h5read(DetectedEventsCTfile, '/beamNo')) + 1; %
beamletNo = double(h5read(DetectedEventsCTfile, '/beamletNo')) + 1; %
energy = h5read(DetectedEventsCTfile, '/energy'); %
eventIds = double(h5read(DetectedEventsCTfile, '/eventIds')) + 1; %
globalTimes = h5read(DetectedEventsCTfile, '/globalTimes'); %

det_dim1 = max(detectorIds1);
det_dim2 = max(detectorIds2);
detectorIds = (detectorIds2-1)*det_dim1+ detectorIds1;

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

clearvars Ind_coin1 Ind_coin2 mask_sameenergy sortInd_sameparticle sortAlldetectorIDInd sortedAlldetectorID timediff badIDbuff1 badIDbuff2

% boundary issue(det id: 1 and 1440)
AlldetectorID2 = (AlleventID-1)*nb_cryst + mod(detectorIds,nb_cryst);
[sortedAlldetectorID, sortAlldetectorIDInd] = sort(AlldetectorID2);
sortInd_sameparticle = find(diff(sortedAlldetectorID)==1);
Ind_coin1 = sortAlldetectorIDInd(sortInd_sameparticle);
Ind_coin2 = sortAlldetectorIDInd(sortInd_sameparticle+1);
mask_sameenergy = (energy(Ind_coin1)-energy(Ind_coin2)==0);
timediff = globalTimes(Ind_coin1)-globalTimes(Ind_coin2);
badIDbuff1 = Ind_coin1(mask_sameenergy & timediff>0);
badIDbuff2 = Ind_coin2(mask_sameenergy & timediff<0);
badID2 = union(badIDbuff1,badIDbuff2);

badID = union(badID1,badID2);
clearvars Ind_coin1 Ind_coin2 mask_sameenergy sortInd_sameparticle sortAlldetectorIDInd sortedAlldetectorID timediff badIDbuff1 badIDbuff2 badID1 badID2

beamletIDs(badID,:) = [];
detectorIds(badID,:) = [];
beamNo(badID,:) = [];
beamletNo(badID,:) = [];
energy(badID,:) = [];
eventIds(badID,:) = [];
globalTimes(badID,:) = [];
Mega(badID,:) = [];
% save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes','-v7.3');

%% fluence map segments
numbeamlets = size(M,2);

x_ = ones(numbeamlets,1);
dose = reshape(M*x_,size(masks{1}.mask));
figure;imshow3D(dose,[])

selectbeamNo = 2;
x_onebeam = zeros(numbeamlets,1);
x_onebeam(cumsumbeamSizes(selectbeamNo)+1:cumsumbeamSizes(selectbeamNo+1)) = 1;
dose_onebeam = reshape(M*x_onebeam,size(masks{1}.mask));
figure;imshow3D(dose_onebeam,[])

% %% Time correction: Adding time of previous events
% doserate = 0.1/60; % (0.1Gy/min)
% time = max(dose_onebeam(:))*numevent/doserate;
% eventrate = time/numevent*1e+09; % ns/event
% 
% numeventbatch = 1e+06;
% numbatches = ceil(numevent/numeventbatch);
% deltatime_event = normrnd(eventrate,eventrate/5,numeventbatch,numbeamlets);
% cumsum_eventtime_batch = cumsum(deltatime_event);
% 
% event_perbatchIDs = mod(eventIds-1, numeventbatch)+1;
% event_batchID = (eventIds - event_perbatchIDs)/numeventbatch + 1;
% event_perbatch_beamletIDs = sub2ind([numeventbatch,numbeamlets], event_perbatchIDs, beamletIDs);
% 
% numbeams = max(beamNo(:));
% batchtime = max(cumsum_eventtime_batch(event_perbatch_beamletIDs));
% beamtime = batchtime*numbatches;
% 
% CorrectedTime = globalTimes + cumsum_eventtime_batch(event_perbatch_beamletIDs)...
%         + event_batchID*batchtime + beamNo*beamtime;
% [sortedtime, sortInd] = sort(CorrectedTime);
% save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'CorrectedTime','sortedtime','sortInd','numevent','eventrate','-append');


%%
Projection_beam = full(sparse(detectorIds,beamNo,1,det_dim1*det_dim2,max(beamNo)));
Projection_beam = reshape(Projection_beam,[det_dim1,det_dim2,size(Projection_beam,2)]);

figure;imshow3D(Projection_beam,[]);


