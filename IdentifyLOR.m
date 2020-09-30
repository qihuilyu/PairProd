function Ind_coin = IdentifyLOR(energy, sortedtime, sortInd, CoincidenceTime, EnergyResolution)

Ind_accept = find(abs(energy-0.511)<0.511*EnergyResolution);

sortInd_coin = find(diff(sortedtime)<CoincidenceTime);
Ind_coin1 = sortInd(sortInd_coin);
Ind_coin2 = sortInd(sortInd_coin+1);

[~, iInd_coin1, ~] = intersect(Ind_coin1, Ind_accept);
[~, iInd_coin2, ~] = intersect(Ind_coin2, Ind_accept);

iInd_coin = intersect(iInd_coin1,iInd_coin2);
Ind_coin1_accept = Ind_coin1(iInd_coin);
Ind_coin2_accept = Ind_coin2(iInd_coin);

Ind_coin = [Ind_coin1_accept Ind_coin2_accept];

