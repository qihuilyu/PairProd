clear
close all
clc
% 
% t = tic;
% while toc(t)<7200*3
%     pause(2)
% end
% !('/media/raid1/qlyu/PairProd/datatest/collect_doescalc_pairprod.sh')

patientName = 'CTphantom_20beams_50m';
projectName = 'PairProd';
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
h5file = fullfile(dosecalcFolder,'PairProd_beamletdose_run2.h5');
maskfile = fullfile(dosecalcFolder,'PairProd_masks.h5');
fmapsfile = fullfile(dosecalcFolder,'PairProd_fmaps.h5');
Anni3Dfile = fullfile(dosecalcFolder,'PairProd_NofPositronAnni3D_run2.h5');
DetectedEventsfile = fullfile(dosecalcFolder,'PairProd_DetectedEvents_run2.h5');


%% masks, fmaps, dose matrix, annihilation matrix
[M,M_Anni,dose_data,masks]=BuildDoseMatrix_PairProd(h5file, maskfile, fmapsfile, Anni3Dfile);  

dose = reshape(full(sum(M,2)),size(masks{1}.mask));
Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));

figure;imshow3D(dose)
figure;imshow3D(Anni3D)

%%

filename = '/media/raid1/qlyu/PairProd/datatest/CTphantom_20beams_50m/dosecalc/test.txt';
% A = readmatrix(filename);
% A = A(:,1);
fileID = fopen(filename,'r');
formatSpec = '%d';
A = fscanf(fileID,formatSpec);
fclose(fileID);

%%

nzbeamlets = find(A);
zbeamlets = find(A==0);

test = full(sum(M,1)*1e+08);
nzInd = find(test);
zInd = find(test==0);


%%

for ii = 1:numel(zbeamlets)
    dose1beam = reshape(full(M(:,zbeamlets(ii))),size(masks{1}.mask));
    tet(ii) = nnz(dose1beam(:));
    
    dose1beam = reshape(full(M_Anni(:,zbeamlets(ii))),size(masks{1}.mask));
    tet2(ii) = nnz(dose1beam(:));
end


%%
%%

for ii = 1:numel(nzbeamlets)
    dose1beam = reshape(full(M(:,nzbeamlets(ii))),size(masks{1}.mask));
    tet22221(ii) = nnz(dose1beam(:));
    
        dose1beam = reshape(full(M_Anni(:,nzbeamlets(ii))),size(masks{1}.mask));
    tet22222(ii) = nnz(dose1beam(:));
end
    
% figure;imshow3D(dose1beam)









