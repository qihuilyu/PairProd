clear
close all
clc

patientName = 'phantom_nanoparticles_50m';
projectName = 'PairProd';
patFolder = fullfile('D:\datatest\PairProd\',patientName);
projectFolder = fullfile(patFolder,projectName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
dosematrixFolder = fullfile(projectFolder,'dosematrix');
resultFolder = fullfile(projectFolder,'results');
mkdir(resultFolder)

load(fullfile(dosematrixFolder,[patientName projectName '_ringdetection.mat']),'detectorIds','beamNo','beamletNo','energy','eventIds','globalTimes','CorrectedTime','sortedtime','sortInd');
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
iso = CenterOfMass*imgres;

beamangles = zeros(numbeams,1);
for BeamNo = 1:numbeams  
    beamangles(BeamNo) = dose_data.beam_metadata(BeamNo).beam_specs.gantry_rot_rad;
end

%% Compute fluence of full FOV
img_fluence_filled = zeros(nx,ny,numbeams);
for BeamNo = 1:numbeams    
    theta = mod(-beamangles(BeamNo) + 3.1436,2*pi);
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
beamlist = 1:numbeams;
for BeamNo = beamlist
    theta = mod(-beamangles(BeamNo) + 3.1436,2*pi);
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
    if(nnz(maskbeam)>1e+04 && FOV<50)
        maskbeam = (sign1.*sign2>0)';
    end
%     figure;imshow([dose1beam/max(dose1beam(:)) maskbeam])
%     C = imfuse(dose1beam/max(dose1beam(:)),maskbeam,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%     figure(49);imshow(C)
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);

    fluence = fluence + img_fluence_filled(:,:,BeamNo).*maskbeam;
    figure(9);imshow(fluence,[])
end

fluence2 = fluence;
fluence2(mask0==0) = 0;
fluence2 = imgaussfilt(fluence2);
figure;imshow(fluence2,[])


Ind = BeamletInd(:,:,beamlist);
Ind = Ind(Ind>0);
xf = zeros(numbeamlets,1);
xf(Ind) = numevents;
doseselectbeam3D = reshape(full(M*xf),size(masks{1}.mask));
doseselectbeam = doseselectbeam3D(:,:,slicenum);
doseselectbeam(mask0==0) = 0;
figure;imshow(doseselectbeam*10,[]);colormap(jet);colorbar; set(gca,'FontSize',30)

Anni1selectbeam3D = reshape(full(M_Anni*xf),size(masks{1}.mask));
Anni1selectbeam = Anni1selectbeam3D(:,:,slicenum);
figure;imshow(Anni1selectbeam,[])

Anni1selectbeam_corrected = Anni1selectbeam./fluence2;
Anni1selectbeam_corrected(mask0==0) = 0;
Anni1selectbeam_corrected(fluence2<max(fluence2(:))*0.2) = 0;
figure(10);imshow(Anni1selectbeam_corrected,[])

%%

ROIInd = [[80.5 115.5 7 8]
        [35.5 84.5 4 7]
        [48.0625 102.3125 4.125 6.25]
        [68.5 108.5 4 6] 
    ];
ROIrow1 = ceil(ROIInd(:,2));
ROIcolumn1 = ceil(ROIInd(:,1));
ROIrow2 = floor(ROIInd(:,2))+floor(ROIInd(:,4));
ROIcolumn2 = floor(ROIInd(:,1))+floor(ROIInd(:,3));

for i = 1
    for j = 1:size(ROIInd,1)
        ImgROI = Anni1selectbeam_corrected(ROIrow1(j):ROIrow2(j),ROIcolumn1(j):ROIcolumn2(j),i);
        if(j==1)
            imginten0 = mean(ImgROI(:));
        end
        ImgROI = ImgROI/imginten0;
        imginten(j) = mean(ImgROI(:));
        imgnoise(j) = std(ImgROI(:));
    end
end

test = (imginten-imginten(1))/imginten(1)*100;
figure;scatter([73,79,83],test(2:end))

xlabel('Atomic No')
ylabel('Increased contrast to water (%)')
set(gca,'FontSize',15)
refline


%% Identify LOR
EnergyResolution = 0.1;
CoincidenceTime = 1;  % ns

Ind_coin_511 = IdentifyLOR_511(energy, CorrectedTime, CoincidenceTime);
Ind_coin_accept = IdentifyLOR(energy, CorrectedTime, CoincidenceTime, EnergyResolution);

TruePositive = length(Ind_coin_511)/length(Ind_coin_accept);
% save(fullfile(dosematrixFolder,[patientName projectName '_detid_pair.mat']),'Ind_coin_511','Ind_coin_accept');

%% Image Reconstruction
TimeResolution = 0.4; % 400 ps
CorrectedTime_TR = CorrectedTime + TimeResolution*randn(size(CorrectedTime));
Ind_coin_accept = IdentifyLOR(energy, CorrectedTime_TR, CoincidenceTime, EnergyResolution);

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
mumapnew = 0*mumap;
mumapnew(ind1+1:end,ind2+1:end) = mumapold(1:end-ind1,1:end-ind2 );
% mumap(1:end-ind1,ind2+1:end) = mumapold(ind1+1:end,1:end-ind2 );
figure;imshow(mumapnew,[])

G = Gtomo2_strip(sg, ig);
% ci = GetACfactor_sino(G, mumap);
li = G * mumapnew;
ci = exp(-li);
img_fbp = em_fbp_QL(sg, ig, sino./ci);

%% Reconstruction-less image generation
TimeResolution = 0.02; % 20 ps
CorrectedTime_TR = CorrectedTime + TimeResolution*randn(size(CorrectedTime));
Ind_coin_accept = IdentifyLOR(energy, CorrectedTime_TR, CoincidenceTime, EnergyResolution);

detid_pair = detectorIds(Ind_coin_accept);
cilist = GetACfactor_list(ci, detid_pair, nb_cryst, R1, distrange, newunidist);
reconparams = struct('nb_cryst',nb_cryst,'R1',R1,'distrange',distrange,...
    'imgres',imgres,'imgsize',[ig.nx, ig.ny]);
deltat = CorrectedTime_TR(Ind_coin_accept(:,1)) - CorrectedTime_TR(Ind_coin_accept(:,2));
[img_direct, img_ci, img_ciN] = recon_TOF_direct(reconparams, detid_pair, deltat, cilist);
% [img_direct, img_ci, img_ciN] = recon_TOF_direct(reconparams, detid_pair, deltat, cilist./cilist);
figure;imshow([Anni2D/max(Anni2D(:)) img_fbp/max(img_fbp(:)) img_direct/max(img_direct(:))],[0.2,1])
Anni2D_corrected = Anni2D./fluence2;
Anni2D_corrected(mask0==0)=0;
img_fbp_corrected = img_fbp./fluence2;
img_fbp_corrected(mask0==0)=0;
img_direct_corrected = img_direct./fluence2;
img_direct_corrected(mask0==0)=0;
figure;imshow([Anni2D_corrected/max(Anni2D_corrected(:)) img_fbp_corrected/max(img_fbp_corrected(:)) img_direct_corrected/max(img_direct_corrected(:))],[0.2,1])

