function [uniang, indang, unidist, inddist, gooddetind] = list2sino(detid_pair, numdet, R1, distrange)

theta = detid_pair/numdet*2*pi;
x_ = R1*sin(theta);
y_ = R1*cos(theta);

p1 = [x_(:,1), y_(:,1)];
p2 = [x_(:,2), y_(:,2)];
pC = (p1+p2)/2;
dist = round(sqrt(sum(pC.^2,2)),4);

p1_p2 = p1 - p2;
xdy = p1_p2(:,1)./p1_p2(:,2);
ang = mod(round(atan(xdy),4),3.1416);

gooddetind = find((~isnan(ang)) & dist<distrange);
ang = ang(gooddetind);
dist = dist(gooddetind);
pC = pC(gooddetind,:);
[uniang, ~, indang] = unique(ang);

xC = pC(:,1);
yC = pC(:,2);
dist(yC<0) = -dist(yC<0);
dist(yC==0) = dist(yC==0).*sign(xC(yC==0));
[unidist, ~, inddist] = unique(dist);

