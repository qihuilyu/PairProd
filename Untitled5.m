%% Compute fluence for limited FOV

beamletwidth = 5;
beamletpix = beamletwidth/imgres;
FOV = size(BeamletLog0,1)*beamletwidth;
beamletlist = 41:numbeamlets;
for BeamletNo = beamletlist
    [xlet,ylet,BeamNo] = ind2sub(size(BeamletInd),find(BeamletInd==BeamletNo));
    theta = mod(-beamangles(BeamNo) + 3.1436,2*pi);
    src = 1000*[cos(theta) sin(theta) 0] + iso;
    
    xf = zeros(numbeamlets,1);
    xf(BeamletNo) = numevents;
    dose1beam3D = reshape(full(M*xf),size(masks{1}.mask));
    dose1beam = dose1beam3D(:,:,slicenum);
    Anni1beam3D = reshape(full(M_Anni*xf),size(masks{1}.mask));
    Anni1beam = Anni1beam3D(:,:,slicenum);
       
    upperlim = (FOV/2-(xlet-1)*beamletwidth)*[-sin(theta) cos(theta) 0] + iso;
    lowerlim = (FOV/2-(xlet)*beamletwidth)*[-sin(theta) cos(theta) 0] + iso;
    
    [y,x] = ndgrid((1:ny)*imgres,(1:nx)*imgres);
    sign1 = sign(y-src(2)-(upperlim(2) - src(2))/(upperlim(1) - src(1))*((x - src(1))));
    sign2 = sign(y-src(2)-(lowerlim(2) - src(2))/(lowerlim(1) - src(1))*((x - src(1))));
    maskbeam = (sign1.*sign2<=0)';
    if(nnz(maskbeam)>1e+04 && FOV<50)
        maskbeam = (sign1.*sign2>0)';
    end
%     figure(100);imshow([dose1beam/max(dose1beam(:)) maskbeam])
    C = imfuse(dose1beam/max(dose1beam(:)),maskbeam,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    figure(49);imshow(C)

%     fluence = fluence + img_fluence_filled(:,:,BeamNo).*maskbeam;
%     figure(9);imshow(fluence,[])
end

% fluence2 = fluence;
% fluence2(mask0==0) = 0;
% fluence2 = imgaussfilt(fluence2);
% figure;imshow(fluence2,[])
% 
% Ind = BeamletInd(:,:,beamletlist);
% Ind = Ind(Ind>0);
% xf = zeros(numbeamlets,1);
% xf(Ind) = numevents;
% doseselectbeam3D = reshape(full(M*xf),size(masks{1}.mask));
% doseselectbeam = doseselectbeam3D(:,:,slicenum);
% doseselectbeam(mask0==0) = 0;
% figure;imshow(doseselectbeam*10,[]);colormap(jet);colorbar; set(gca,'FontSize',30)
% 
% Anni1selectbeam3D = reshape(full(M_Anni*xf),size(masks{1}.mask));
% Anni1selectbeam = Anni1selectbeam3D(:,:,slicenum);
% figure;imshow(Anni1selectbeam,[])
% 
% Anni1selectbeam_corrected = Anni1selectbeam./fluence2;
% Anni1selectbeam_corrected(mask0==0) = 0;
% Anni1selectbeam_corrected(fluence2<max(fluence2(:))*0.2) = 0;
% figure(10);imshow(Anni1selectbeam_corrected,[])
