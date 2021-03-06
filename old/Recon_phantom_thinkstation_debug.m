clear
close all
clc

patientName = 'phantom_tumorwithAuCa_1m_10MV';
projectName = 'PairProd';
patFolder = fullfile('D:\datatest\PairProd\',patientName);
projectFolder = fullfile(patFolder,projectName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
dosematrixFolder = fullfile(projectFolder,'dosematrix');
resultFolder = fullfile(projectFolder,'result');
mkdir(resultFolder)

load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'energy','detectorIds','CorrectedTime','eventIds');
load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');

paramsFolder = fullfile(projectFolder,'params');
ParamsNum = 0;
load(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
InfoNum = 0;
load(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');

% slicenum = 88;
slicenum = 50;
x_CT = img(:,:,end+1-slicenum)-1000;
[mumap,densmap,Ind] = lookup_materials_bulk_density(x_CT);

% figure;imshow(mumap,[])
% figure;imshow(densmap,[])
% figure;imshow(Ind,[])
% figure;imshow(x_CT,[])

%% Identify LOR
EnergyResolution = 0.1;
CoincidenceTime = 2;  % ns 

Ind_coin_511 = IdentifyLOR_511(energy, CorrectedTime, CoincidenceTime);
Ind_coin_accept = IdentifyLOR(energy, CorrectedTime, CoincidenceTime, EnergyResolution);

TruePositive = length(Ind_coin_511)/length(Ind_coin_accept);
% save(fullfile(dosematrixFolder,[patientName projectName '_detid_pair.mat']),'Ind_coin_511','Ind_coin_accept');

%% Image Reconstruction
R1 = 1200;
distrange = 300;
imgsize = size(img);
nb_cryst = max(detectorIds);

detid_pair = detectorIds(Ind_coin_accept);
[sino, dr, newunidist, sinobuff, unidist] = rebin_PET2(detid_pair, nb_cryst, R1, distrange);

ig = image_geom('nx', size(img,1), 'ny', size(img,2), 'fov', size(img,1)*imgres);
sg = sino_geom('par', 'nb', size(sino,1), 'na', size(sino,2), 'dr', dr);
img_fbp_nocorrect = em_fbp_QL(sg, ig, sino);

Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
Anni2D = Anni3D(:,:,slicenum);
figure;imshow(Anni2D,[])

dose3D = reshape(full(sum(M,2)),size(masks{1}.mask));
dose2D = dose3D(:,:,slicenum);
figure;imshow(dose2D,[])

Anni2D = Anni3D(:,:,slicenum);
Anni2Dold = Anni2D;
Anni2D = 0*Anni2D;
ind1 = 20; ind2 = 20;
Anni2D(ind1+1:end,ind2+1:end) = Anni2Dold(1:end-ind1,1:end-ind2 );
% Anni2D(1:end+ind1,ind2+1:end) = Anni2Dold(-ind1+1:end,1:end-ind2 );
% figure;imshow([Anni2D/max(Anni2D(:)),img_fbp_nocorrect/max(img_fbp_nocorrect(:))],[])
C = imfuse(Anni2D/max(Anni2D(:)),img_fbp_nocorrect/max(img_fbp_nocorrect(:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(C)
% li = G * mumap;
% ci = exp(-li);
% img_fbp = em_fbp_QL(sg, ig, sino./ci);
% figure;imshow([img_fbp_gtsino/max(img_fbp_gtsino(:)) img_fbp/max(img_fbp(:))],[])

mumapold = mumap;
mumap = 0*mumap;
mumap(ind1+1:end,ind2+1:end) = mumapold(1:end-ind1,1:end-ind2 );
% mumap(1:end-ind1,ind2+1:end) = mumapold(ind1+1:end,1:end-ind2 );
figure;imshow(mumap,[])



G = Gtomo2_strip(sg, ig);
% ci = GetACfactor_sino(G, mumap);
li = G * mumap;
ci = exp(-li);
figure;imshow([sino/max(sino(:)); ci/max(ci(:))])
figure;imshow([sino/max(sino(:)); li/max(li(:))])
img_fbp = em_fbp_QL(sg, ig, sino./ci);
img_fbp_nocorrect = em_fbp_QL(sg, ig, sino);
figure; imshow([sino/max(sino(:))],[0.1,0.8])
figure; imshow([ci/max(ci(:))],[0,0.01])

ForBack.applyFP = @(x) G*x(:);
ForBack.applyBP = @(x) G'*x;
gamma = 2000;
mu = 1e-05;
[x_TV, Maincost_TV] = IterRecon_PairProd_TV_FISTA (ForBack, sino./ci, gamma, mu, [ig.nx, ig.ny]);

% proj = (G * img_fbp) .* ci;
% figure;imshow([proj; Sino],[])


sino_FP = reshape(ForBack.applyFP(Anni2D),size(sino));
img_fbp_gtsino = em_fbp_QL(sg, ig, sino_FP);
figure;imshow(sino_FP,[])
sino_correct = sino./ci;

test = sino_correct/sum(sino_correct(:))-sino_FP/sum(sino_FP(:));
figure;imshow([test/sum(abs(test(:))); sino_correct/sum(sino_correct(:)); sino_FP/sum(sino_FP(:))],[])
figure;imshow([sino_correct/sum(sino_correct(:)); sino_FP/sum(sino_FP(:))],[])
figure;imshow([sinobuff/sum(sinobuff(:)); sino/sum(sino(:)); sino_FP/sum(sino_FP(:))],[])
figure;imshow([img_fbp_gtsino/max(img_fbp_gtsino(:)) img_fbp/max(img_fbp(:)) img_fbp_nocorrect/max(img_fbp_nocorrect(:))],[])
figure;imshow([sino_correct/sum(sino_correct(:))-sino_FP/sum(sino_FP(:))],[])
figure;imshow([img_fbp_nocorrect/max(img_fbp_nocorrect(:))],[])
figure;imshow([img_fbp/max(img_fbp(:))],[])

figure;imshow([img_fbp_gtsino/max(img_fbp_gtsino(:))-img_fbp_nocorrect/max(img_fbp_nocorrect(:))],[])

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
figure;imshow(img_ci/max(img_ci(:)),[0,0.005])
% figure;imshow([Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:)) img_fbp/max(img_fbp(:)) ],[])
 figure;imshow([[Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:))]; [img_direct/max(img_direct(:)) x_TV/max(x_TV(:))]],[])
 figure;imshow([Anni2D/max(Anni2D(:)) img_fbp_nocorrect/max(img_fbp_nocorrect(:)) img_direct2/max(img_direct2(:))],[])

%% TOF reconstruction
count = 1;
for TR = [0.2,0.6,1,1.5 2]
    [sino3D, dr, sinobuff3D] = rebin_PET_TOF(reconparams, detid_pair, deltat, TR);
    sino3DAll{count} = sino3D;
    sinobuff3DAll{count} = sinobuff3D;
    [nd,na] = size(ci);
    ci3D = repmat(reshape(ci,[nd,1,na]),[1,nd,1]);
    sino3DAll{count} = flip(flip(sino3D,1),2);
    img_toffbp{count} = em_toffbp_QL(sg, ig, sino3DAll{count});
    img_tofbp{count} = em_tof_backprojection_QL(sg, ig, sino3DAll{count});
   
    figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_toffbp{count}/max(img_toffbp{count}(:)) img_tofbp{count}/max(img_tofbp{count}(:)) img_direct/max(img_direct(:))],[])
    title(['TR: ' num2str(TR) 's'])
    saveas(gcf,fullfile(resultFolder,['imgrecon_TR' num2str(TR) 'corrected.png']))
    
    count = count+1;
end



figure;
subplot(3,1,1); imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_direct/max(img_direct(:))],[])
subplot(3,1,2); imshow([img_toffbp{1}/max(img_toffbp{1}(:)) img_toffbp{2}/max(img_toffbp{2}(:)) img_toffbp{3}/max(img_toffbp{3}(:)) img_toffbp{4}/max(img_toffbp{4}(:))],[])
subplot(3,1,3); imshow([img_tofbp{1}/max(img_tofbp{1}(:)) img_tofbp{2}/max(img_tofbp{2}(:)) img_tofbp{3}/max(img_tofbp{3}(:)) img_tofbp{4}/max(img_tofbp{4}(:))],[])
saveas(gcf,fullfile(resultFolder,['Compare_TR.png']))

