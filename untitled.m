clear
close all
clc

patientName = 'phantom_10cm_mingap_10k';
projectName = 'PairProd';
patFolder = fullfile('/media/raid1/qlyu/PairProd/datatest',patientName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
h5file = fullfile(dosecalcFolder,'PairProd_beamletdose.h5');
maskfile = fullfile(dosecalcFolder,'PairProd_masks.h5');
fmapsfile = fullfile(dosecalcFolder,'PairProd_fmaps.h5');
Anni3Dfile = fullfile(dosecalcFolder,'PairProd_NofPositronAnni3D.h5');
DetectedEventsfile = fullfile(dosecalcFolder,'PairProd_DetectedEvents.h5');

%% masks, fmaps, dose matrix, annihilation matrix
[M,M_Anni,dose_data,masks]=BuildDoseMatrix_PairProd(h5file, maskfile, fmapsfile, Anni3Dfile);  

dose = reshape(full(sum(M,2)),size(masks{1}.mask));
Anni3D = reshape(full(sum(M_Anni,2)),size(masks{1}.mask));

figure;imshow3D(dose)
figure;imshow3D(Anni3D)

pdose = 25;
[StructureInfo, params] = InitIMRTparams_DLMCforRyan(M,dose_data,masks,pdose,[1,2,3]);


%% 
filename = '/media/raid1/qlyu/PairProd/testcode/dosecalc/deepmc_debug_data/deepmc_debug_data/mask.npy';
data = readNPY(filename);
figure;imshow3D(double(data),[])
