load('D:\datatest\PairProd\phantom_nanoparticles_2mmbeamlet_5m\PairProd\results\Recon_pairprod.mat')
load('D:\datatest\PairProd\phantom_nanoparticles_360beam_460m_merged\CTsim\results\Recon_CT.mat')
i = 1;
ImgInfo(i).Img_raw = Anni2D.*mask0;
ImgInfo(i).Img_corrected = Anni2D_corrected.*mask0;
ImgInfo(i).ImgNorMax = ImgInfo(i).Img_raw/max(ImgInfo(i).Img_raw(:));
ImgInfo(i).Method = 'PPI ground-truth';
i = i+1;
ImgInfo(i).Img_raw = img_fbp.*mask0;
ImgInfo(i).Img_corrected = img_fbp_corrected.*mask0;
ImgInfo(i).ImgNorMax = ImgInfo(i).Img_raw/max(ImgInfo(i).Img_raw(:));
ImgInfo(i).Method = 'PPI FBP';
i = i+1;
ImgInfo(i).Img_raw = img_beampath.*mask0;
ImgInfo(i).Img_corrected = img_beampath_corrected.*mask0;
ImgInfo(i).ImgNorMax = ImgInfo(i).Img_raw/max(ImgInfo(i).Img_raw(:));
ImgInfo(i).Method = 'PPI beam-path';
i = i+1;
ImgInfo(i).Img_raw = img_direct.*mask0;
ImgInfo(i).Img_corrected = img_direct_corrected.*mask0;
ImgInfo(i).ImgNorMax = ImgInfo(i).Img_raw/max(ImgInfo(i).Img_raw(:));
ImgInfo(i).Method = 'PPI reconstruction-less';
i = i+1;
ImgInfo(i).Img_raw = CT_FBP(:,:,end/2).*mask0;
ImgInfo(i).Img_corrected = CT_FBP(:,:,end/2).*mask0;
ImgInfo(i).ImgNorMax = ImgInfo(i).Img_raw/max(ImgInfo(i).Img_raw(:));
ImgInfo(i).Method = 'CT FBP';

figure;imshow([ImgInfo.ImgNorMax])


%%

ROIInd = [[99.1875 82.5625 8.00000000000001 7.87499999999999]
        [102.5 63.5 5 5]
        [89.5 46.5 5 5]
        [68.5 39.5 5 5]
        [48.5 46.5 4 5]
        [35.5 62.5 4 7]
        [35.5 84.5 4 7]
        [48.0625 102.3125 4.125 6.25]
        [68.5 108.5 4 6] 
    ]; 

ROIrow1 = ceil(ROIInd(:,2));
ROIcolumn1 = ceil(ROIInd(:,1));
ROIrow2 = floor(ROIInd(:,2))+floor(ROIInd(:,4));
ROIcolumn2 = floor(ROIInd(:,1))+floor(ROIInd(:,3));

for i = 1:numel(ImgInfo)
    Img = ImgInfo(i).Img;
    for j = 1:size(ROIInd,1)
        ImgROI = Img(ROIrow1(j):ROIrow2(j),ROIcolumn1(j):ROIcolumn2(j));
        
        if(j==1)
            imginten0 = mean(ImgROI(:));
        end
        ImgROI = ImgROI/imginten0;
        imginten(i,j) = mean(ImgROI(:));
        imgnoise(i,j) = std(ImgROI(:));
    end
end


%%
test = (imginten-imginten(:,1))./repmat(imginten(:,1),[1,size(imginten,2)])*100;
AtomicNum = [53,56,64,70,73,79,83];
MarkerSize = 50;
LineWidth = 2;
figure; 
% yyaxis right;
% scatter(AtomicNum,test(4,3:end),MarkerSize,'s','filled');  hold on;
% ylabel('Increased contrast to water (%), CT')
% legend({'CT'})
scatter(AtomicNum,test(1,3:end),MarkerSize,'o','filled'); hold on;
scatter(AtomicNum,test(2,3:end),MarkerSize,'+','LineWidth',LineWidth);  hold on;
scatter(AtomicNum,test(3,3:end),MarkerSize,'d','LineWidth',LineWidth);  hold on;
ylabel('Increased contrast to water (%)')
legend({'Ground truth','Reconstruction less','FBP'})


% icolor = 'r';
% figure;scatter(AtomicNum,test(1,3:end),[],icolor);hline=refline; hline.Color = icolor;  hold on;
% icolor = 'y';
% scatter(AtomicNum,test(2,3:end),[],icolor);hline2=refline; hline2.Color = icolor;
% icolor = 'b';
% scatter(AtomicNum,test(3,3:end),[],icolor);hline3=refline; hline3.Color = icolor;  hold on;

xlabel('Atomic No')
set(gca,'FontSize',15)


figure; 
scatter(AtomicNum,test(4,3:end),MarkerSize,'s','filled');  hold on;
ylabel('Increased contrast to water (%)')
legend({'CT'})
xlabel('Atomic No')
set(gca,'FontSize',15)


%%


figure;imagesc([CT_Dose(:,:,end/2) PP_Dose(:,:,end/2)*10]); colorbar; colormap(jet)



