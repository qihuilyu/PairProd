clear
close all
clc

patientName = 'GBMHY';
projectName = 'PairProd';
patFolder = fullfile('D:\datatest\PairProd\',patientName);
projectFolder = fullfile(patFolder,projectName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
dosematrixFolder = fullfile(projectFolder,'dosematrix');
resultFolder = fullfile(projectFolder,'result');
mkdir(resultFolder)

load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes','CorrectedTime','sortedtime','sortInd');
load(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','M_Anni','dose_data','masks');
load(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');


%%


load('D:\datatest\PairProd\GBMHY\PairProd\results\GBMHY PairProd Info7 eta0 result.mat')
xPolish = result.xPolish;
StructureInfo = result.StructureInfo;

slicenum = 88;
x_CT = img(:,:,end+1-slicenum)-1000;
[mumap,densmap,Ind] = lookup_materials_bulk_density(x_CT);

% figure;imshow(mumap,[])
% figure;imshow(densmap,[])
% figure;imshow(Ind,[])
% figure;imshow(x_CT,[])

Anni3D = reshape(full(M_Anni*xPolish),size(StructureInfo(1).Mask));
Anni2D = Anni3D(:,:,slicenum);

dose3D = reshape(full(M*xPolish),size(StructureInfo(1).Mask));
dose2D = dose3D(:,:,slicenum);
figure;imshow(dose2D,[])

Anni2D = Anni3D(:,:,slicenum);
Anni2Dold = Anni2D;
Anni2D = 0*Anni2D;
ind1 = 20; ind2 = 20;
Anni2D(ind1+1:end,ind2+1:end) = Anni2Dold(1:end-ind1,1:end-ind2 );

mumapold = mumap;
mumap = 0*mumap;
mumap(ind1+1:end,ind2+1:end) = mumapold(1:end-ind1,1:end-ind2 );
% mumap(1:end-ind1,ind2+1:end) = mumapold(ind1+1:end,1:end-ind2 );
figure;imshow(mumap,[])


%% Identify LOR
EnergyResolution = 0.1;
CoincidenceTime = 3;  % ns 

Ind_coin_511 = IdentifyLOR_511(energy, sortedtime, sortInd, CoincidenceTime);
Ind_coin_accept = IdentifyLOR(energy, sortedtime, sortInd, CoincidenceTime, EnergyResolution);

TruePositive = length(Ind_coin_511)/length(Ind_coin_accept);
% save(fullfile(dosematrixFolder,[patientName projectName '_detid_pair.mat']),'Ind_coin_511','Ind_coin_accept');

%% Image Reconstruction
R1 = 1200;
distrange = 300;
imgsize = size(img);
nb_cryst = max(detectorIds);

detid_pair = detectorIds(Ind_coin_accept);
beamSizes = squeeze(sum(sum(params.BeamletLog0,1),2));
cumsumbeamSizes = cumsum([0; beamSizes]);
beamNoshift = cumsumbeamSizes(beamNo);
beamletIDs = double(beamletNo) + beamNoshift;
beamletID_select = beamletIDs(Ind_coin_accept(:,1));
optimizedweights = xPolish(beamletID_select);
[sino, dr, newunidist, sinobuff, unidist] = rebin_PET2(detid_pair, nb_cryst, R1, distrange, optimizedweights);

ig = image_geom('nx', size(img,1), 'ny', size(img,2), 'fov', size(img,1)*imgres);
sg = sino_geom('par', 'nb', size(sino,1), 'na', size(sino,2), 'dr', dr);
G = Gtomo2_strip(sg, ig);
% ci = GetACfactor_sino(G, mumap);
li = G * mumap;
ci = exp(-li);
img_fbp = em_fbp_QL(sg, ig, sino./ci);
img_fbp_nocorrect = em_fbp_QL(sg, ig, sino);

ForBack.applyFP = @(x) G*x(:);
ForBack.applyBP = @(x) G'*x;
gamma = 2000;
mu = 1e-05;
[x_TV, Maincost_TV] = IterRecon_PairProd_TV_FISTA (ForBack, sino./ci, gamma, mu, [ig.nx, ig.ny]);


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
[img_direct, img_ci, img_ciN] = recon_TOF_direct(reconparams, detid_pair, deltat, cilist./optimizedweights);
img_direct2 = recon_TOF_direct(reconparams, detid_pair, deltat, 1./optimizedweights);
% profile report
% figure; imshow([img_direct/max(img_direct(:)) img_direct2/max(img_direct2(:)) img_ci/max(img_ci(:)) img_ciN/max(img_ciN(:))],[])
figure;imshow(img_ci/max(img_ci(:)),[0,0.005])
% figure;imshow([Anni2D/max(Anni2D(:)) img_direct/max(img_direct(:)) img_fbp/max(img_fbp(:)) ],[])
 figure;imshow([[Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:))]; [img_direct/max(img_direct(:)) x_TV/max(x_TV(:))]],[])
 figure;imshow([Anni2D/max(Anni2D(:)) img_fbp_nocorrect/max(img_fbp_nocorrect(:)) img_direct2/max(img_direct2(:))],[])

%%
img_fbp_nocorrect(img_fbp_nocorrect<0) = 0;
img_fbp_nocorrect = img_fbp_nocorrect/max(img_fbp_nocorrect(:));

planName = 'img_fbp_nocorrect';
xoff = -5;
yoff = 5;
addDoseToGui_Move_QL(repmat(img_fbp_nocorrect,[1,1,size(Anni3D,3)]),[planName],xoff,yoff)


planName = 'img_direct';
xoff = -5;
yoff = 5;
addDoseToGui_Move_QL(repmat(img_direct,[1,1,size(Anni3D,3)]),[planName],xoff,yoff)


Anni3D = reshape(full(M_Anni*xPolish),size(StructureInfo(1).Mask));
planName = 'Anni3D';
xoff = -0.2;
yoff = 0.2;
addDoseToGui_Move_QL(Anni3D,[planName 'anni'],xoff,yoff)

