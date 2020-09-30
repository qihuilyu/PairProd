clear
close all
clc

patientName = 'phantom_10cm';
projectName = 'PairProd';
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
projectFolder = fullfile(patFolder,projectName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
dosematrixFolder = fullfile(projectFolder,'dosematrix');

load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']));
load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');

paramsFolder = fullfile(projectFolder,'params');
ParamsNum = 0;
load(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
InfoNum = 0;
load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');

mask = StructureInfo(2).Mask(:,:,56);
mumap = double(mask);
mumap(mask==1) = 0.0096;
mumap(mask==0) = 0.0001;

%% Identify LOR
EnergyResolution = 0.1;
CoincidenceTime = 2;  % ns 

Ind_coin_511 = IdentifyLOR_511(energy, sortedtime, sortInd, CoincidenceTime);
Ind_coin_accept = IdentifyLOR(energy, sortedtime, sortInd, CoincidenceTime, EnergyResolution);

TruePositive = length(Ind_coin_511)/length(Ind_coin_accept);

%% Image Reconstruction
R1 = 1200;
distrange = 300;
imgsize = size(img);
nb_cryst = max(detectorIds);

detid_pair = detectorIds(Ind_coin_accept);
[sino, dr, sinobuff] = rebin_PET(detid_pair, nb_cryst, R1, distrange);
figure;imshow(sino,[])
% figure;imshow(sino./ci,[])
% figure;imshow(ci,[])


ig = image_geom('nx', imgsize(1), 'ny', imgsize(2), 'fov', imgsize(1)*imgres);
sg = sino_geom('par', 'nb', size(sino,1), 'na', nb_cryst, 'dr', dr);
figure; sg.plot([ig])	
G = Gtomo2_strip(sg, ig);
ci = GetACfactor_sino(G, mumap);

img_fbp = em_fbp_QL(sg, ig, flip(sino./ci,1));
% figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_direct2/max(img_direct2(:))],[])


ForBack.applyFP = @(x) G*x;
ForBack.applyBP = @(x) G'*x;
gamma = 500;
mu = 1e-05;
[x_TV, Maincost_TV] = IterRecon_PairProd_TV_FISTA (ForBack, sino(:)./ci(:), gamma, mu, [ig.nx, ig.ny]);

img_fbp2 = em_fbp(sg, ig, sino./ci);
img_fbp3 = em_fbp_QL(sg, ig, sino./ci);
img_fbp2 = em_fbp_QL(sg, ig, sino./ci);

% proj = (G * img_fbp) .* ci;
% figure;imshow([proj; Sino],[])

Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
Anni2D = Anni3D(:,:,ceil(end/2));
figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:))],[])
figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_fbp2/max(img_fbp2(:)) img_fbp3/max(img_fbp3(:))],[])


%% Reconstruction-less image generation
cilist = GetACfactor_list(ci, detid_pair, nb_cryst, R1, distrange);
reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
    'imgres',imgres,'imgsize',size(img_fbp));
profile on
img_direct = recon_TOF_direct(CorrectedTime, detectorIds, Ind_coin_accept, cilist./cilist, reconparams);
profile report
img_direct2 = recon_TOF_direct(CorrectedTime, detectorIds, Ind_coin_accept, cilist, reconparams);
figure; imshow(img_direct,[])
figure;imshow([Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:)) img_fbp/max(img_fbp(:)) ],[])
figure;imshow([Anni2D/max(Anni2D(:)) img_fbp2/max(img_fbp2(:)) img_direct2/max(img_direct2(:))],[])

detectorIds_360 = floor(detectorIds/4);
reconparams_360 = struct('nb_cryst',max(detectorIds_360),'R1',R1,'distrange',distrange,...
    'imgres',imgres,'imgsize',size(img_fbp));
profile on
img_direct3 = recon_TOF_direct(CorrectedTime, detectorIds_360, Ind_coin_accept, cilist./cilist, reconparams_360);
profile report
figure;imshow([Anni2D/max(Anni2D(:)) img_direct3/max(img_direct3(:)) img_direct/max(img_direct(:))],[])

%% TOF reconstruction
[sino3D, dr, sinobuff3D] = rebin_PET_TOF(CorrectedTime, detectorIds, Ind_coin_accept, reconparams);
[nd,na] = size(ci);
ci3D = repmat(reshape(ci,[nd,1,na]),[1,nd,1]);
sinoflipped = flip(flip(sino3D./ci3D,1),2);
[img_toffbp] = em_toffbp_QL(sg, ig, sinoflipped);
 
 figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_toffbp/max(img_toffbp(:)) img_direct/max(img_direct(:))],[])






