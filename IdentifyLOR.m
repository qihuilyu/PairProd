function Ind_coin = IdentifyLOR(energy, sortedtime, sortInd, CoincidenceTime, EnergyResolution)

Ind_accept = find(abs(energy-0.511)<0.511*EnergyResolution);

difftime = diff(sortedtime);
difftime2 = difftime(1:end-1) + difftime(2:end);
sortInd_coin = find(difftime<CoincidenceTime);
sortInd_coin2 = find(difftime2<CoincidenceTime);
sortInd_multiplecoin = union(union(sortInd_coin2,sortInd_coin2+1),sortInd_coin2+2);
sortInd_coin_cleaned = setdiff(sortInd_coin,sortInd_multiplecoin);

Ind_coin1 = sortInd(sortInd_coin_cleaned);
Ind_coin2 = sortInd(sortInd_coin_cleaned+1);

[~, iInd_coin1, ~] = intersect(Ind_coin1, Ind_accept);
[~, iInd_coin2, ~] = intersect(Ind_coin2, Ind_accept);

iInd_coin = intersect(iInd_coin1,iInd_coin2);
Ind_coin1_accept = Ind_coin1(iInd_coin);
Ind_coin2_accept = Ind_coin2(iInd_coin);

Ind_coin = [Ind_coin1_accept Ind_coin2_accept];

