spec60kVp = [
15.0e3  3.69755e1
 15.5e3  7.60284e1
 16.0e3  1.44202e2
 16.5e3  2.50543e2
 17.0e3  4.14450e2
 17.5e3  6.39255e2
 18.0e3  9.45618e2
 18.5e3  1.34175e3
 19.0e3  1.82652e3
 19.5e3  2.41447e3
 20.0e3  3.09233e3
 20.5e3  3.84197e3
 21.0e3  4.71803e3
 21.5e3  5.62049e3
 22.0e3  6.62026e3
 22.5e3  7.63616e3
 23.0e3  8.71456e3
 23.5e3  9.74511e3
 24.0e3  1.08834e4
 24.5e3  1.19128e4
 25.0e3  1.29030e4
 25.5e3  1.39615e4
 26.0e3  1.49513e4
 26.5e3  1.58474e4
 27.0e3  1.67786e4
 27.5e3  1.75912e4
 28.0e3  1.84313e4
 28.5e3  1.91184e4
 29.0e3  1.98201e4
 29.5e3  2.03437e4
 30.0e3  2.08710e4
 30.5e3  2.14021e4
 31.0e3  2.17308e4
 31.5e3  2.21404e4
 32.0e3  2.24149e4
 32.5e3  2.26680e4
 33.0e3  2.28511e4
 33.5e3  2.29787e4
 34.0e3  2.30580e4
 34.5e3  2.31126e4
 35.0e3  2.30857e4
 35.5e3  2.30500e4
 36.0e3  2.29651e4
 36.5e3  2.28233e4
 37.0e3  2.26729e4
 37.5e3  2.24744e4
 38.0e3  2.22677e4
 38.5e3  2.20139e4
 39.0e3  2.17295e4
 39.5e3  2.14377e4
 40.0e3  2.11162e4
 40.5e3  2.07658e4
 41.0e3  2.03945e4
 41.5e3  1.99959e4
 42.0e3  1.95918e4
 42.5e3  1.91819e4
 43.0e3  1.87333e4
 43.5e3  1.82994e4
 44.0e3  1.78223e4
 44.5e3  1.73535e4
 45.0e3  1.68564e4
 45.5e3  1.63673e4
 46.0e3  1.58573e4
 46.5e3  1.53279e4
 47.0e3  1.48167e4
 47.5e3  1.42818e4
 48.0e3  1.37346e4
 48.5e3  1.31995e4
 49.0e3  1.26480e4
 49.5e3  1.20818e4
 50.0e3  1.15261e4
 50.5e3  1.09583e4
 51.0e3  1.03896e4
 51.5e3  9.81948e3
 52.0e3  9.23999e3
 52.5e3  8.66863e3
 53.0e3  8.08867e3
 53.5e3  7.50874e3
 54.0e3  6.92891e3
 54.5e3  6.34320e3
 55.0e3  5.76421e3
 55.5e3  5.18055e3
 56.0e3  4.59827e3
 56.5e3  4.01745e3
 57.0e3  3.43817e3
 57.5e3  2.86051e3
 58.0e3  2.28457e3
 58.5e3  1.70884e3
 59.0e3  1.13721e3
 59.5e3  5.67110e2
 60.0e3  -4.72342e1];

figure;scatter(spec60kVp(:,1),spec60kVp(:,2),'LineWidth',2);
title('60kVp spectrum')

nn = 256;
spec60kVp0 = [spec60kVp; [ones(nn-size(spec60kVp,1),1)*spec60kVp(end,1) zeros(nn-size(spec60kVp,1),1)]];

spec60kVpDS = [];
spec60kVpDS(:,2) = spec60kVp0(1:2:end,2)+spec60kVp0(2:2:end,2);
spec60kVpDS(:,1) = (spec60kVp0(1:2:end,1)+spec60kVp0(2:2:end,1))/2;
figure;scatter(spec60kVpDS(:,1),spec60kVpDS(:,2),'LineWidth',2);
title('60kVp spectrum fist Downsample')

spec60kVp0 = spec60kVpDS; spec60kVpDS = [];
spec60kVpDS(:,2) = spec60kVp0(1:2:end,2)+spec60kVp0(2:2:end,2);
spec60kVpDS(:,1) = (spec60kVp0(1:2:end,1)+spec60kVp0(2:2:end,1))/2;
figure;scatter(spec60kVpDS(:,1),spec60kVpDS(:,2),'LineWidth',2);
title('60kVp spectrum second Downsample')

% spec60kVp0 = spec60kVpDS; spec60kVpDS = [];
% spec60kVpDS(:,2) = spec60kVp0(1:2:end,2)+spec60kVp0(2:2:end,2);
% spec60kVpDS(:,1) = (spec60kVp0(1:2:end,1)+spec60kVp0(2:2:end,1))/2;
% figure;scatter(spec60kVpDS(:,1),spec60kVpDS(:,2),'LineWidth',2);
% title('60kVp spectrum third Downsample')

spec60kVpmask = spec60kVpDS(:,2)>0;
spec60kVpfinal = spec60kVpDS(spec60kVpmask==1,:);

spec60kVpfinal(:,1) = spec60kVpfinal(:,1)/1e+06;
spec60kVpfinal(:,2) = spec60kVpfinal(:,2)/sum(spec60kVpfinal(:,2));
figure;scatter(spec60kVpfinal(:,1),spec60kVpfinal(:,2),'LineWidth',2);
title('60kVp spectrum final')

