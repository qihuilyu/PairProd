clear
close all
clc

patientName = 'GBMHY2';
projectName = 'PairProd';
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
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
    if(~strcmp(dcinfo.Modality,'RTSTRUCT'))
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
badID1 = Ind_coin1(mask_sameenergy & timediff>0);
badID2 = Ind_coin2(mask_sameenergy & timediff<0);
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

deltatime_event = normrnd(eventrate,eventrate/5,numevent,numbeamlets);
cumsum_eventtime = cumsum(deltatime_event);
event_beamletIDs = sub2ind([numevent,numbeamlets], eventIds, beamletIDs);

numbeams = max(beamNo(:));
beamtime = max(cumsum_eventtime(event_beamletIDs));
deltatime_beam = beamtime*(1:numbeams)'; % 1 ms

CorrectedTime = globalTimes + cumsum_eventtime(event_beamletIDs) + deltatime_beam(beamNo);
[sortedtime, sortInd] = sort(CorrectedTime);
save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'CorrectedTime','sortedtime','sortInd','-append');
