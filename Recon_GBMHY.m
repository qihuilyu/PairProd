clear
close all
clc

patientName = 'GBMHY2';
projectName = 'PairProd';
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
projectFolder = fullfile(patFolder,projectName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
dosematrixFolder = fullfile(projectFolder,'dosematrix');

load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']));
load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');

%% Identify LOR
CoincidenceTime = 2;  % ns  
Ind_coin_511 = IdentifyLOR_511(energy, sortedtime, sortInd, CoincidenceTime);

EnergyResolution = 0.1;
Ind_coin_accept = IdentifyLOR(energy, sortedtime, sortInd, CoincidenceTime, EnergyResolution);
TruePositive = length(Ind_coin_511)/length(Ind_coin_accept);

%% Image Reconstruction
R1 = 120;
distrange = 30;
imgsize = size(img);
nb_cryst = max(detectorIds);

detid_pair = detectorIds(Ind_coin_accept);
Sino = rebin_PET(detid_pair, nb_cryst, R1, distrange);
figure;imshow(Sino,[])

ntheta = size(Sino,1);
dr = 2*distrange/ntheta;

ig = image_geom('nx', imgsize(1), 'ny', imgsize(2), 'fov', imgsize(1)*imgres);
sg = sino_geom('par', 'nb', ntheta, 'na', nb_cryst, 'dr', dr);
img_fbp = em_fbp(sg, ig, Sino);
figure;imshow(img_fbp/max(img_fbp(:)),[])

Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
Anni2D = Anni3D(:,:,ceil(end/2));
figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:))],[])


%% Reconstruction-less image generation
reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
    'imgres',imgres,'imgsize',size(img_fbp));

img_direct = PETrecon_TOF_direct(CorrectedTime, detectorIds, Ind_coin_accept, reconparams);
figure; imshow(img_direct,[])
figure;imshow([img_fbp/max(img_fbp(:)) Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:))],[])




