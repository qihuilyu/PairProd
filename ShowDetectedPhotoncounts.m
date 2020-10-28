Ind_511 = find(abs(energy-0.511)<0.00001);
EnergyResolution = 0.1;
Ind_accept = find(abs(energy-0.511)<0.511*EnergyResolution);

beamdet_511 = full(sparse(beamNo(Ind_511),detectorIds(Ind_511),1));
figure; set(gcf,'pos',[2715   148    1480    1001])
for ii = 1:2:size(beamdet_511,1)
    plot(beamdet_511(ii,:),'LineWidth',2); hold on
end
% legend({'beam 1','beam 5','beam 9','beam 13','beam 17'})
legend({'beam 1','beam 3','beam 5','beam 7'})
xlabel('Detector ID')
ylabel('Detected photon counts')
title('Detected 511keV photon counts')
set(gca,'FontSize',20)

beamdet_accept = full(sparse(beamNo(Ind_accept),detectorIds(Ind_accept),1));
figure; set(gcf,'pos',[2715   148    1480    1001])
for ii = 1:2:size(beamdet_accept,1)
    plot(beamdet_accept(ii,:),'LineWidth',2); hold on
end
% legend({'beam 1','beam 5','beam 9','beam 13','beam 17'})
legend({'beam 1','beam 3','beam 5','beam 7'})
xlabel('Detector ID')
ylabel('Detected photon counts')
title('Detected photon counts within 10% energy resolution')
set(gca,'FontSize',20)


thresh = median(beamdet_accept(ii,:))*3;
[row,col] = ind2sub(size(beamdet_accept),find(beamdet_accept>thresh));

