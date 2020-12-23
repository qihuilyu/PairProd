clear
close all
clc

% t = tic;
% while toc(t)<7200
%     pause(2)
% end
% !('/media/raid1/qlyu/PairProd/datatest/collect_doescalc_pairprod.sh')

patientName = 'GBMHY_final_100m';
projectName = 'PairProd';
patFolder = fullfile('/media/raid0/qlyu/PairProd/datatest',patientName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
maskfile = fullfile(dosecalcFolder,'PairProd_masks.h5');
fmapsfile = fullfile(dosecalcFolder,'PairProd_fmaps.h5');
numeventsfile = fullfile(dosecalcFolder,[patientName '_nruns_10.txt']);
taglist = {'run01','run02','run03','run04','run05','run06','run07','run08','run09','run10'};

projectFolder = fullfile(patFolder,'PairProd');
paramsFolder = fullfile(projectFolder,'params');
mkdir(paramsFolder)
dosematrixFolder = fullfile(projectFolder,'dosematrix');
mkdir(dosematrixFolder)

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

%% numevents
fileID = fopen(numeventsfile,'r');
formatSpec = '%d';
numeventsvec = fscanf(fileID,formatSpec);
fclose(fileID);

%% Multiple runs
for itag = 1:numel(taglist)
    tag = taglist{itag};
    
    ParamsNum = 0;
    load(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
    InfoNum = 0;
    load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');
    
    load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes_' tag '.mat']),'M','M_Anni','dose_data','masks');
        
    load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection_' tag '.mat']),'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes');
    
    numevent = max(eventIds);
    numevent = numevent + 99 - mod(numevent-1,100);
    beamSizes = squeeze(sum(sum(params.BeamletLog0,1),2));
    cumsumbeamSizes = cumsum([0; beamSizes]);
    beamNoshift = cumsumbeamSizes(beamNo);
    beamletIDs = double(beamletNo) + beamNoshift;
    
    %% fluence map segments
    numbeamlets = size(M,2);
    
    x_ = numeventsvec;
    dose = reshape(M*x_,size(masks{1}.mask));
    figure;imshow3D(dose,[])
    
    
    %% Time correction: Adding time of previous events
    load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'eventrate');
    eventrate = eventrate*10;
    
    numeventbatch = 1e+06;
    numbatches = ceil(numevent/numeventbatch);
    deltatime_event = normrnd(eventrate,eventrate,numeventbatch,numbeamlets);
    cumsum_eventtime_batch = cumsum(deltatime_event);
    
    event_perbatchIDs = mod(eventIds-1, numeventbatch)+1;
    event_batchID = (eventIds - event_perbatchIDs)/numeventbatch + 1;
    event_perbatch_beamletIDs = sub2ind([numeventbatch,numbeamlets], event_perbatchIDs, beamletIDs);
    
    numbeams = max(beamNo(:));
    batchtime = max(cumsum_eventtime_batch(event_perbatch_beamletIDs));
    beamtime = batchtime*numbatches;
    
    CorrectedTime = globalTimes + cumsum_eventtime_batch(event_perbatch_beamletIDs)...
        + event_batchID*batchtime + beamNo*beamtime;
    [sortedtime, sortInd] = sort(CorrectedTime);
    save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection_' tag '_detectorefficiency0p1.mat']),...
        'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes','beamletIDs','CorrectedTime','sortedtime','sortInd','numeventsvec');
    
end