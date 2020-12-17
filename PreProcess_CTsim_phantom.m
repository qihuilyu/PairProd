clear
close all
clc
% 
% t = tic;
% while toc(t)<7200*3
%     pause(2)
% end
% !('/media/raid1/qlyu/PairProd/datatest/collect_doescalc_pairprod.sh')

patientName = 'phantom_nanoparticles_360beam_50m_thinslice5mm_CTsimNEW_run3';
projectName = 'CTsim';
Folder = '/media/raid1/qlyu/PairProd/datatest';
patFolder = fullfile(Folder,patientName);
dosecalcFolder = fullfile(patFolder,'dosecalc');
h5file = fullfile(dosecalcFolder,[projectName '_beamletdose.h5']);
maskfile = fullfile(dosecalcFolder,[projectName '_masks.h5']);
fmapsfile = fullfile(dosecalcFolder,[projectName '_fmaps.h5']);
DetectedEventsCTfile = fullfile(dosecalcFolder,[projectName '_CTprojection.mat']);

%% masks, fmaps, dose matrix, annihilation matrix
[M,dose_data,masks]=BuildDoseMatrix_CTsim(h5file, maskfile, fmapsfile);

dose = reshape(full(sum(M,2)),size(masks{1}.mask));
figure;imshow3D(dose)

pdose = 25;
[StructureInfo, params] = InitIMRTparams_DLMCforRyan(M,dose_data,masks,pdose,[1,2,0]);

projectFolder = fullfile(patFolder,projectName);
paramsFolder = fullfile(projectFolder,'params');
mkdir(paramsFolder)

ParamsNum = 0;
save(fullfile(paramsFolder,['params' num2str(ParamsNum) '.mat']),'params');
InfoNum = 0;
save(fullfile(paramsFolder,['StructureInfo' num2str(InfoNum) '.mat']),'StructureInfo');

dosematrixFolder = fullfile(projectFolder,'dosematrix');
mkdir(dosematrixFolder)
save(fullfile(dosematrixFolder,[patientName projectName '_M_HighRes.mat']),'M','dose_data','masks','-v7.3');

%% Dicom
DicomPath = fullfile(dosecalcFolder,'ctdata');
baseFileNames = dir([DicomPath '/*.dcm']);
[sortedFile,index] = sort_nat({baseFileNames.name});

img = [];
for ii = 1:numel(sortedFile)
    dcinfo = dicominfo(fullfile(DicomPath,sortedFile{ii}));
    if(~strcmpi(dcinfo.Modality,'RTstruct'))
        imgres = dcinfo.SliceThickness;
        img(:,:,ii) = dicomread(fullfile(DicomPath,sortedFile{ii}));
    end
end
figure;imshow3D(img,[0,2000])
save(fullfile(dosematrixFolder,[patientName projectName '_dicomimg.mat']),'img','imgres');

%%
load(fullfile(Folder,'phantom_nanoparticles_360beam_200m_thinslice5mm_CTsimNEW_blankfield','dosecalc','CTsim_CTprojection.mat'),'CTprojection');
Projection0 = CTprojection(:,:,1)/4;

load(DetectedEventsCTfile,'CTprojection');
LI = log(permute(repmat(Projection0,[1,1,size(CTprojection,3)]),[2,1,3])./permute(CTprojection,[2,1,3]));

LI(isinf(LI)) = 0;
LI(isnan(LI)) = 0;
figure;imshow3D(LI,[]);
save(fullfile(dosematrixFolder,'LineIntegrals.mat'),'CTprojection','LI');


%%


	cg = ct_geom('fan', ...
	'ns', 250, ... % detector channels
	'nt', 500, ... % detector rows
	'na', 360, ... % angular samples
	'offset_s', 0, ... % quarter-detector offset
	'dsd', 1000, ...
	'dod', 333.3, ...
	'dfs', inf, ... % arc
	'ds', 2, ... % detector pitch
	'dt', 2, ... % detector row spacing for 0.625mm slices, 2009-12-06
	'pitch',0,...
    'orbit_start',90);
	ig = image_geom('nx', size(img,1), 'ny', size(img,2), 'nz', 100, 'fov', size(img,1)*imgres);
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
    li_hat = fdk_filter(LI(126:end-125,:,1:end), 'ramp', cg.dsd, cg.dfs, cg.ds);
	
    
    args = {flip(li_hat,1), cg, ig, 'ia_skip', 1}; % increase 1 for faster debugging
	CT_FBP = cbct_back(args{:}, 'use_mex', 1, 'back_call', @jf_mex);
    figure;imshow3D(CT_FBP,[0,0.15])
    
 resultsFolder = fullfile(projectFolder,'results');
 mkdir(resultsFolder)
 save(fullfile(resultsFolder,'Recon_CT.mat'),'CT_FBP')
   
% 	back2 = cbct_back(args{:}, 'use_mex', 1, 'back_call', @fdk_mex);
% 	max_percent_diff(back1, back2)
% 	xfdk = feldkamp(cg, ig, li_hat, ...
% 		'extrapolate_t', ceil(1.3 * cg.nt/2)); % todo: compute carefully
% 
% 	pr [ig.nx ig.ny ig.nz ig.dx ig.dy ig.dz]
% 	pr [cg.ns cg.nt cg.na cg.ds cg.dt]
% 	printm('image Mbytes %d', ig.nx*ig.ny*ig.nz*4 / 1024^2)
% 	printm('proj Mbytes %d', cg.ns*cg.nt*cg.na*4 / 1024^2)
% 
% 	printm 'ellipse proj' % somewhat realistic phantom object
% 	ell = [ ...
% 		[30 10 10	150 150 280	0 0 1000]; % 30cm diam
% 		[80 10 10	50 50 30	0 0 300]; % bone-like
% 		[-10 -40 75	40 40 40	0 0 600];
% 		[-10 80 -20	30 30 30	0 0 600];
% 	];
% %	xtrue = ellipsoid_im(ig, ell); im(xtrue); return
% 
% 	proj = ellipsoid_proj(cg, ell);
% 	proj = fdk_filter(proj, 'ramp', cg.dsd, cg.dfs, cg.ds);
% 	if 0 % zero outer edges of projection for testing
% 		proj([1 cg.ns], :, :) = 0;
% 		proj(:, [1 cg.nt], :) = 0;
% 	end



