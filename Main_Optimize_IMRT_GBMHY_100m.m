% clearvars -except M
% close all
% clc
% g = gpuDevice(1)
% reset(g)

patientName = 'GBMHY_100m';
folder = 'D:\datatest\PairProd\';
patFolder = fullfile(folder, patientName);
bz2Folder = fullfile(patFolder,'CERR_bz2');
OutputFileName = fullfile(bz2Folder,[patientName '.mat']);
% CERR('CERRSLICEVIEWER')
% sliceCallBack_QL('OPENNEWPLANC', OutputFileName);

projectName = 'PairProd';
projectFolder = fullfile(patFolder,[projectName]);
dosematrixFolder = fullfile(projectFolder,'dosematrix');
paramsFolder = fullfile(projectFolder,'params');
resultsFolder = fullfile(projectFolder,'results');
mkdir(resultsFolder)

load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']))

InfoNum = 7;
load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']));
% save(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');
% StructureInfo(2).Mask = imfill(StructureInfo(2).Mask,'holes');

paramsNum = 0;
load(fullfile(paramsFolder,['params' num2str(paramsNum) '.mat']));

%% Downsample and prepare matrix
DS = 1;
[A,Weights] = CreateA_DLMCforRyan(M, StructureInfo);
ATrans = A';

[Dx,Dy] = CreateDxDyFMO(params.BeamletLog0);
D = [Dx;Dy];

%% beam selection
seed = 2; % Seed for the random number generator.
rng(seed) % Set random number generator seed so we can reproduce experiments
tic
params.stepSize = 1e-04;
params.beamWeight = 5000;
params.maxIter = 8000;
params.showTrigger = 100;

PTV = StructureInfo(1).Mask;
BODY = StructureInfo(2).Mask;
finalBeams = 1:numel(params.beamSizes);

%% Polish step
paramsPolish = params;
eta = 0;
paramsPolish.eta = eta;
paramsPolish.maxIter = 600;
tic
[xPolish,costsDF_polish,costs_polish] = polish_BOO_IMRT_cpu(finalBeams,A,D,Weights,paramsPolish);
timePolish = toc;
figure;semilogy(costs_polish)
% xPolish(Ind) = 0;
xPolish = gather(xPolish);
fsize = size(params.BeamletLog0);
xfull = zeros(fsize(1),fsize(2),16);
xfull(params.BeamletLog0==1) = xPolish;
xfull4D = reshape(xfull,[size(xfull,1),size(xfull,2),4,4]);
xfull2D = reshape(permute(xfull4D,[1,3,2,4]),[size(xfull,1)*4,size(xfull,2)*4]);
figure;imshow(xfull2D,[]); colormap(jet); colorbar;
figureFolder = fullfile(patFolder, 'figures');
mkdir(figureFolder)
saveas(gcf,fullfile(figureFolder,['fmap ' patientName ' ' projectName ' Info' num2str(InfoNum)...
    ' params' num2str(paramsNum) '.png']))

%% Visualize and save results
dose = M*xPolish; dose = reshape(dose,size(PTV)); dose(BODY==0&PTV==0)=0;
planName = [patientName ' ' projectName ' Info' num2str(InfoNum) ' eta' num2str(eta)];
xoff = -0.2;
yoff = 0.2;
% addDoseToGui_DLMCforRyan(BODY,[planName],xoff,yoff,zoff);
% addDoseToGui_dvo(dose,[planName]);
addDoseToGui_Move_QL(dose,[planName],xoff,yoff)

N_anni = M_Anni*xPolish;
N_anni = reshape(N_anni,size(PTV)); N_anni(BODY==0&PTV==0)=0;
planName = [patientName ' ' projectName ' Info' num2str(InfoNum) ' eta' num2str(eta)];
xoff = -0.2;
yoff = 0.2;
% addDoseToGui_DLMCforRyan(BODY,[planName],xoff,yoff,zoff);
% addDoseToGui_dvo(dose,[planName]);
addDoseToGui_Move_QL(N_anni/max(N_anni(:))*max(dose(:)),[planName 'anni'],xoff,yoff)

load(fullfile(patFolder,[patientName '_DoseInfo.mat']),'DoseInfo');
PlanIndex = length(DoseInfo)+1;
DoseInfo(PlanIndex).Data = dose;
DoseInfo(PlanIndex).Name = planName;
DoseInfo(PlanIndex).DFCost = costsDF_polish(end);
DoseInfo(PlanIndex).Cost = costs_polish(end);
DoseInfo(PlanIndex).Date = datestr(datetime);
DoseInfo(PlanIndex).xPolish = xPolish;
save(fullfile(patFolder,[patientName '_DoseInfo.mat']),'DoseInfo','-v7.3');

strNum = [1,3:numel(StructureInfo)-2];
numBins = 200;
scale = plotDVH_QL(DoseInfo([end]), strNum, StructureInfo, numBins, 0);

result = struct('patientName',patientName,'projectName',projectName,...
    'dose',dose,'finalBeams',finalBeams,'xPolish',xPolish,...
    'timePolish',timePolish,'costs_polish',costs_polish,...
    'params',params,'StructureInfo',StructureInfo,...
    'paramsPolish',paramsPolish,'planName',planName);
save(fullfile(resultsFolder,[planName ' result.mat']),'result')

for ii = 1:length(StructureInfo)
    StructureInfo(ii).meandose = mean(dose(StructureInfo(ii).Mask==1));
end


