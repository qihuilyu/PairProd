clear
close all
clc

patientName = 'phantom_LimitedROI_50m';
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

slicenum = 50;
x_CT = img(:,:,end+1-slicenum)-1000;
[mumap,densmap,Ind] = lookup_materials_bulk_density(x_CT);

% figure;imshow(mumap,[])
% figure;imshow(densmap,[])
% figure;imshow(Ind,[])
% figure;imshow(x_CT,[])

%%
numevents = max(eventIds);
numevents = numevents + 100 - mod(numevents,100);

Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
Anni2D = Anni3D(:,:,slicenum);
figure;imshow(Anni2D,[])

mumap_10MV = 0*x_CT;
mumap_10MV(x_CT>=0) = 0.004942;
mumap_10MV(x_CT<0) = 0;

PTV = StructureInfo(1).Mask;
BODY = StructureInfo(2).Mask ==1 | StructureInfo(1).Mask ==1;
mask0 = BODY(:,:,slicenum);
figure;imshow(mask0,[])

BeamletLog0 = params.BeamletLog0;
BeamletInd = BeamletLog0;
numbeamlets = nnz(BeamletLog0);
BeamletInd(BeamletLog0==1) = 1:numbeamlets;
numbeams = size(BeamletLog0,3);
[nx,ny,~] = size(img);
[CenterOfMass] = GetPTV_COM(PTV);
iso = CenterOfMass*imgres - [0,1,0];

%% Compute fluence of full FOV
img_fluence_filled = zeros(nx,ny,numbeams);
for BeamNo = 1:numbeams    
    theta = 1.033*pi - BeamNo/numbeams*2*pi + 0.003;
    src = 1000*[cos(theta) sin(theta) 0] + iso;
    tardir = [-sin(theta) cos(theta) 0];
    
    img_fluence= ComputeFluence(mumap_10MV, src, tardir, iso, imgres);
    img_fluence(img_fluence==0) = NaN;
    img_fluence_filled(:,:,BeamNo) = inpaint_nans(img_fluence);
    figure(8);imshow(img_fluence_filled(:,:,BeamNo),[])
end


%% Compute fluence for limited FOV

fluence = 0*mumap_10MV;
FOV = size(BeamletLog0,1)*15/imgres-1;
beamlist = 1:60
for BeamNo = beamlist
    theta = 1.033*pi - BeamNo/numbeams*2*pi + 0.003;
    src = 1000*[cos(theta) sin(theta) 0] + iso;

    Ind = BeamletInd(:,:,BeamNo);
    Ind = Ind(Ind>0);
    
    xf = zeros(numbeamlets,1);
    xf(Ind) = numevents;
    dose1beam3D = reshape(full(M*xf),size(masks{1}.mask));
    dose1beam = dose1beam3D(:,:,slicenum);
    Anni1beam3D = reshape(full(M_Anni*xf),size(masks{1}.mask));
    Anni1beam = Anni1beam3D(:,:,slicenum);
       
    upperlim = -FOV*[-sin(theta) cos(theta) 0] + iso;
    lowerlim = FOV*[-sin(theta) cos(theta) 0] + iso;
    
    [y,x] = ndgrid((1:ny)*imgres,(1:nx)*imgres);
    sign1 = sign(y-src(2)-(upperlim(2) - src(2))/(upperlim(1) - src(1))*((x - src(1))));
    sign2 = sign(y-src(2)-(lowerlim(2) - src(2))/(lowerlim(1) - src(1))*((x - src(1))));
    maskbeam = (sign1.*sign2<0)';
    if(nnz(maskbeam)>1e+04)
        maskbeam = (sign1.*sign2>0)';
    end
    % figure;imshow([dose1beam/max(dose1beam(:)) maskbeam'])
%     C = imfuse(dose1beam/max(dose1beam(:)),maskbeam,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%     figure(49);imshow(C)
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);

    fluence = fluence + img_fluence_filled(:,:,BeamNo).*maskbeam;
    figure(9);imshow(fluence,[])
end

fluence2 = fluence;
fluence2(mask0==0) = 0;
fluence2 = imgaussfilt(fluence2);


Ind = BeamletInd(:,:,beamlist);
Ind = Ind(Ind>0);
xf = zeros(numbeamlets,1);
xf(Ind) = numevents;
doseselectbeam3D = reshape(full(M*xf),size(masks{1}.mask));
doseselectbeam = doseselectbeam3D(:,:,slicenum);
figure;imshow(doseselectbeam,[])

Anni1selectbeam3D = reshape(full(M_Anni*xf),size(masks{1}.mask));
Anni1selectbeam = Anni1selectbeam3D(:,:,slicenum);
figure;imshow(Anni1selectbeam,[])

Anni1selectbeam_corrected = Anni1selectbeam./fluence2;
Anni1selectbeam_corrected(mask0==0) = 0;
Anni1selectbeam_corrected(fluence2<max(fluence2(:))*0.3) = 0;
figure(10);imshow(Anni1selectbeam_corrected,[])

%%

ROIInd = [
        [35.5 62.5 4 7]
        [35.5 84.5 4 7]
        [48.0625 102.3125 4.125 6.25]
    ];
ROIrow1 = ceil(ROIInd(:,2));
ROIcolumn1 = ceil(ROIInd(:,1));
ROIrow2 = floor(ROIInd(:,2))+floor(ROIInd(:,4));
ROIcolumn2 = floor(ROIInd(:,1))+floor(ROIInd(:,3));

for i = 1
    for j = 1:size(ROIInd,1)
        ImgROI = Anni1selectbeam_corrected(ROIrow1(j):ROIrow2(j),ROIcolumn1(j):ROIcolumn2(j),i);
        imginten(j) = mean(ImgROI(:));
        imgnoise(j) = std(ImgROI(:));
    end
end

test = (imginten-imginten(1))/imginten(1)*100;
figure;scatter([53,56,64,70,73,79,83],test(3:end))

xlabel('Atomic No')
ylabel('Increased contrast to water (%)')
set(gca,'FontSize',15)
refline



%% Compute fluence for limited FOV

fluence = 0*mumap_10MV;
FOV = size(BeamletLog0,1)*15/imgres-1;
for BeamNo = 1:numbeams    
    theta = 1.033*pi - BeamNo/numbeams*2*pi + 0.003;
    src = 1000*[cos(theta) sin(theta) 0] + iso;

    Ind = BeamletInd(:,:,BeamNo);
    Ind = Ind(Ind>0);
    
    xf = zeros(numbeamlets,1);
    xf(Ind) = numevents;
    dose1beam3D = reshape(full(M*xf),size(masks{1}.mask));
    dose1beam = dose1beam3D(:,:,slicenum);
    Anni1beam3D = reshape(full(M_Anni*xf),size(masks{1}.mask));
    Anni1beam = Anni1beam3D(:,:,slicenum);
       
    upperlim = -FOV*[-sin(theta) cos(theta) 0] + iso;
    lowerlim = FOV*[-sin(theta) cos(theta) 0] + iso;
    
    [y,x] = ndgrid((1:ny)*imgres,(1:nx)*imgres);
    sign1 = sign(y-src(2)-(upperlim(2) - src(2))/(upperlim(1) - src(1))*((x - src(1))));
    sign2 = sign(y-src(2)-(lowerlim(2) - src(2))/(lowerlim(1) - src(1))*((x - src(1))));
    maskbeam = (sign1.*sign2<0)';
    if(nnz(maskbeam)>1e+04)
        maskbeam = (sign1.*sign2>0)';
    end
    % figure;imshow([dose1beam/max(dose1beam(:)) maskbeam'])
%     C = imfuse(dose1beam/max(dose1beam(:)),maskbeam,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%     figure(49);imshow(C)
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);

    fluence = fluence + img_fluence_filled(:,:,BeamNo).*maskbeam;
%     figure(9);imshow(fluence,[])
end


fluence2 = fluence;
fluence2(mask0==0) = 0;
figure(10);imshow(fluence2,[])

fluence2 = imgaussfilt(fluence2);
figure(10);imshow(fluence2,[])
figure(10);imshow(1./fluence2,[0,0.1])


Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));
Anni2D = Anni3D(:,:,slicenum);
figure;imshow(Anni2D,[])
Anni2D_corrected = Anni2D./fluence2*50000000;
Anni2D_corrected(mask0==0) = 0;
figure(10);imshow(Anni2D_corrected,[102,137])


figure;imshow(Anni2D/max(Anni2D(:)),[])
figure;imshow(fluence2/max(fluence2(:)),[])


C = imfuse(Anni2D/max(Anni2D(:)),fluence/max(fluence(:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure(49);imshow(C)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);






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
figure;imshow(img_fbp_nocorrect,[])

dose3D = reshape(full(sum(M,2)),size(masks{1}.mask));
dose2D = dose3D(:,:,slicenum);
figure;imshow(dose2D,[])

ind1 = 1; ind2 = 0;
Anni2D = Anni3D(:,:,slicenum);
Anni2Dold = Anni2D;
Anni2D = 0*Anni2D;
Anni2D(ind1+1:end,ind2+1:end) = Anni2Dold(1:end-ind1,1:end-ind2 );
% Anni2D(1:end+ind1,ind2+1:end) = Anni2Dold(-ind1+1:end,1:end-ind2 );
% figure;imshow([Anni2D/max(Anni2D(:)),img_fbp_nocorrect/max(img_fbp_nocorrect(:))],[])
C = imfuse(Anni2D/max(Anni2D(:)),img_fbp_nocorrect/max(img_fbp_nocorrect(:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure; imshow(C)

mumapold = mumap;
mumap = 0*mumap;
mumap(ind1+1:end,ind2+1:end) = mumapold(1:end-ind1,1:end-ind2 );
% mumap(1:end-ind1,ind2+1:end) = mumapold(ind1+1:end,1:end-ind2 );
figure;imshow(mumap,[])

G = Gtomo2_strip(sg, ig);
% ci = GetACfactor_sino(G, mumap);
li = G * mumap;
ci = exp(-li);
img_fbp = em_fbp_QL(sg, ig, sino./ci);

%% Reconstruction-less image generation
cilist = GetACfactor_list(ci, detid_pair, nb_cryst, R1, distrange, newunidist);
reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
    'imgres',imgres,'imgsize',[ig.nx, ig.ny]);
deltat = CorrectedTime(Ind_coin_accept(:,1)) - CorrectedTime(Ind_coin_accept(:,2));
[img_direct, img_ci, img_ciN] = recon_TOF_direct(reconparams, detid_pair, deltat, cilist);
figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_direct/max(img_direct(:))],[])

