
filename = '/media/raid1/qlyu/PairProd/dbdata/mcgeometry/5f9a56e8f415784a8a010acf/mcgeo.txt';
A = readmatrix(filename);

% density = reshape(A(:,1),[105,125,65]);
% figure;imshow3D(density,[])
% 
% Ind = reshape(A(:,2),[105,125,65]);
% figure;imshow3D(Ind,[])
% 

ImageSize = [108 128 68];

density = reshape(A(:,1),ImageSize);
figure;imshow3D(density,[])

Ind = reshape(A(:,2),ImageSize);
figure;imshow3D(Ind,[])

