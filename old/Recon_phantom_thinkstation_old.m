clear
close all
clc

patientName = 'phantom_polygon_10cm_10m';
projectName = 'PairProd';
patFolder = fullfile('D:\datatest\PairProd\',patientName);
projectFolder = fullfile(patFolder,projectName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
dosematrixFolder = fullfile(projectFolder,'dosematrix');
resultFolder = fullfile(projectFolder,'result');
mkdir(resultFolder)

load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']));
load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');

paramsFolder = fullfile(projectFolder,'params');
ParamsNum = 0;
load(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
InfoNum = 0;
load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');

hulist =  [-800,   -154,   -18,     197,      1000];  
%          lung  adipose   water   muscle     bone
maclist = [0.0095 0.0096  0.0096   0.0095    0.0089];
rholist = [0.3     0.92      1      1.05      1.6  ];
mulist = maclist.*rholist;

x_CT = img(:,:,ceil(end/2)) - 1024;
mumap = interp1(hulist,mulist,x_CT(:),'linear');
mumap(x_CT>max(hulist)) = max(mulist);
mumap(x_CT<min(hulist)) = min(mulist);
mumap = reshape(mumap,size(x_CT));

figure;imshow(mumap,[])

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
% detectorIds = rem(detectorIds + 180, nb_cryst) + 1;

detid_pair = detectorIds(Ind_coin_accept);
% [sino, dr, sinobuff] = rebin_PET_half(detid_pair, nb_cryst, R1, distrange);
[sino, dr, newunidist, sinobuff, unidist] = rebin_PET2(detid_pair, nb_cryst, R1, distrange);

ig = image_geom('nx', size(img,1), 'ny', size(img,2), 'fov', size(img,1)*imgres);
sg = sino_geom('par', 'nb', size(sino,1), 'na', size(sino,2), 'dr', dr);
G = Gtomo2_strip(sg, ig);
ci = GetACfactor_sino(G, mumap);
figure;imshow([(flip(sino,1))/max(sino(:)); ci/max(ci(:))])
img_fbp = em_fbp_QL(sg, ig, sino./ci);
figure; imshow([sino/max(sino(:))],[0.1,0.8])
figure; imshow([ci/max(ci(:))],[0,0.01])
figure;imshow(proj/max(proj(:)),[0,0.05])


ForBack.applyFP = @(x) G*x(:);
% ForBack.applyBP = @(x) G'*x;
% gamma = 50;
% mu = 1e-05;
% [x_TV, Maincost_TV] = IterRecon_PairProd_TV_FISTA (ForBack, flip(sino,1), gamma, mu, [ig.nx, ig.ny]);

% proj = (G * img_fbp) .* ci;
% figure;imshow([proj/max(proj(:)); sino/max(sino(:))],[])

Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
Anni2D = Anni3D(:,:,ceil(end/2));
% figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:))],[])
% figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_fbp2/max(img_fbp2(:)) img_fbp3/max(img_fbp3(:))],[])

sino_FP = reshape(ForBack.applyFP(Anni2D),size(sino));
img_fbp_gtsino = em_fbp_QL(sg, ig, sino_FP);
figure;imshow(sino_FP,[])
sino_correct = sino./ci;
figure;imshow([sino_correct/sum(sino_correct(:)); sino_FP/sum(sino_FP(:))],[])
figure;imshow([sinobuff/sum(sinobuff(:)); sino/sum(sino(:)); sino_FP/sum(sino_FP(:))],[])
figure;imshow([img_fbp_gtsino/max(img_fbp_gtsino(:)) img_fbp/max(img_fbp(:))],[])
figure;imshow([sino/sum(sino(:))-sino_FP/sum(sino_FP(:))],[])


%% Reconstruction-less image generation
cilist = GetACfactor_list(ci, detid_pair, nb_cryst, R1, distrange, newunidist);
reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
    'imgres',imgres,'imgsize',[ig.nx, ig.ny]);
% profile on
deltat = CorrectedTime(Ind_coin_accept(:,1)) - CorrectedTime(Ind_coin_accept(:,2));
[img_direct, img_ci, img_ciN] = recon_TOF_direct(reconparams, detid_pair, deltat, cilist);
img_direct2 = recon_TOF_direct(reconparams, detid_pair, deltat, cilist./cilist);
% profile report
% figure; imshow([img_direct/max(img_direct(:)) img_direct2/max(img_direct2(:)) img_ci/max(img_ci(:)) img_ciN/max(img_ciN(:))],[])
% figure;imshow([Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:)) img_fbp/max(img_fbp(:)) ],[])

detectorIds_360 = floor(detectorIds/4);
reconparams_360 = struct('nb_cryst',max(detectorIds_360),'R1',R1,'distrange',distrange,...
    'imgres',12,'imgsize',[ig.nx, ig.ny]);
profile on
img_direct_360 = recon_TOF_direct(CorrectedTime, detectorIds_360, Ind_coin_accept, cilist, reconparams_360);
profile report
figure;imshow([Anni2D/max(Anni2D(:)) img_direct_360/max(img_direct_360(:)) img_direct/max(img_direct(:))],[])


detid_pair_360 = detectorIds_360(Ind_coin_accept);
[sino_360, dr, sinobuff] = rebin_PET(detid_pair_360, max(detectorIds_360), R1, distrange);
sg_360 = sino_geom('par', 'nb', size(sino_360,1), 'na', size(sino_360,2), 'dr', dr);
img_fbp_360 = em_fbp_QL(sg_360, ig, (flip(sino_360,1)));
figure;imshow([img_fbp_360/max(img_fbp_360(:)) img_direct_360/max(img_direct_360(:))])

%% TOF reconstruction
count = 1;
for TR = [0.2,0.6,1,1.5 2]
    [sino3D, dr, sinobuff3D] = rebin_PET_TOF(CorrectedTime, detectorIds, Ind_coin_accept, TR, reconparams);
    sino3DAll{count} = sino3D;
    sinobuff3DAll{count} = sinobuff3D;
    [nd,na] = size(ci);
    ci3D = repmat(reshape(ci,[nd,1,na]),[1,nd,1]);
    sino3DAll{count} = flip(flip(sino3D,1),2);
    img_toffbp{count} = em_toffbp_QL(sg, ig, sino3DAll{count});
    img_tofbp{count} = em_tof_backprojection_QL(sg, ig, sino3DAll{count});
   
    figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_toffbp{count}/max(img_toffbp{count}(:)) img_tofbp{count}/max(img_tofbp{count}(:)) img_direct/max(img_direct(:))],[])
    title(['TR: ' num2str(TR) 's'])
    saveas(gcf,fullfile(resultFolder,['imgrecon_TR' num2str(TR) '.png']))
    
    count = count+1;
end



figure;
subplot(3,1,1); imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_direct/max(img_direct(:))],[])
subplot(3,1,2); imshow([img_toffbp{1}/max(img_toffbp{1}(:)) img_toffbp{2}/max(img_toffbp{2}(:)) img_toffbp{3}/max(img_toffbp{3}(:)) img_toffbp{4}/max(img_toffbp{4}(:))],[])
subplot(3,1,3); imshow([img_tofbp{1}/max(img_tofbp{1}(:)) img_tofbp{2}/max(img_tofbp{2}(:)) img_tofbp{3}/max(img_tofbp{3}(:)) img_tofbp{4}/max(img_tofbp{4}(:))],[])
saveas(gcf,fullfile(resultFolder,['Compare_TR.png']))

