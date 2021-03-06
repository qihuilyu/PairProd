clear
close all
clc

projectName = 'PairProd';
patientNameList = {'phantom_nanoparticles_LimitedROI','phantom_nanoparticles_LimitedROI_run2'};

M0 = 0;
M_Anni0 = 0;
detectorIds0 = [];
beamNo0 = [];
beamletNo0 = [];
energy0 = [];
eventIds0 = [];
globalTimes0 = [];
numevent0 = 0;


for ii = 1:numel(patientNameList)
    patientName = patientNameList{ii};
    patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
    projectFolder = fullfile(patFolder,projectName);
    paramsFolder = fullfile(projectFolder,'params');
    dosematrixFolder = fullfile(projectFolder,'dosematrix');
    
    InfoNum = 0;
    load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');
    load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');
    ParamsNum = 0;
    load(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
    
    load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
    
    load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),...
        'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes',...
        'numevent');
    numevent = round(numevent/100)*100;
    M0 = (M0*numevent0 + M*numevent)/(numevent + numevent0);
    M_Anni0 = (M_Anni0*numevent0 + M_Anni*numevent)/(numevent + numevent0);
    detectorIds0 = [detectorIds0;detectorIds];
    beamNo0 = [beamNo0;beamNo];
    beamletNo0 = [beamletNo0;beamletNo];
    energy0 = [energy0;energy];
    eventIds = eventIds + numevent0;
    eventIds0 = [eventIds0;eventIds];
    globalTimes0 = [globalTimes0;globalTimes];
    numevent0 = numevent0 + numevent;
    
    
end


M = M0;
M_Anni = M_Anni0;
detectorIds = detectorIds0;
beamNo = beamNo0;
beamletNo = beamletNo0;
energy = energy0;
eventIds = eventIds0;
globalTimes = globalTimes0;
numevent = numevent0;


%% Save files
patientName = [patientNameList{1},'_merged'];
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
projectFolder = fullfile(patFolder,projectName);
paramsFolder = fullfile(projectFolder,'params');
dosematrixFolder = fullfile(projectFolder,'dosematrix');

mkdir(dosematrixFolder)
mkdir(paramsFolder)

InfoNum = 0;
save(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');
save(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');
ParamsNum = 0;
save(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');

dosematrixFolder = fullfile(projectFolder,'dosematrix');
save(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks','-v7.3');

save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),...
    'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes',...
    'numevent');

beamSizes = squeeze(sum(sum(params.BeamletLog0,1),2));
cumsumbeamSizes = cumsum([0; beamSizes]);
beamNoshift = cumsumbeamSizes(beamNo);
beamletIDs = double(beamletNo) + beamNoshift;

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

numeventbatch = 1e+06;
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
save(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'CorrectedTime','sortedtime','sortInd','numevent','-append');







