spec10MV = [
0.15 0.0239
0.20 0.0505
0.30 0.0659
0.40 0.0489
0.50 0.0620
0.60 0.0673
0.80 0.0814
1.00 0.0742
1.25 0.0622
1.50 0.0823
2.00 0.1160
3.00 0.0966
4.00 0.0597
5.00 0.0405
6.00 0.0376
8.00 0.0254
10.00 0.0055];

figure;plot(spec10MV(:,1),spec10MV(:,2),'LineWidth',2);
title('10MV spectrum')
xlabel('Energy (MV)')
ylabel('Percentage')