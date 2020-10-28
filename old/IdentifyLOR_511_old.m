function Ind_coin = IdentifyLOR_511_old(energy, sortedtime, sortInd, CoincidenceTime)

Ind_511 = find(abs(energy-0.511)<0.0001);

difftime = diff(sortedtime);
difftime2 = difftime(1:end-1) + difftime(2:end);
sortInd_coin = find(difftime<CoincidenceTime);
sortInd_coin2 = find(difftime2<CoincidenceTime);
sortInd_multiplecoin = union(union(sortInd_coin2,sortInd_coin2+1),sortInd_coin2+2);
sortInd_coin_cleaned = setdiff(sortInd_coin,sortInd_multiplecoin);

Ind_coin1 = sortInd(sortInd_coin_cleaned);
Ind_coin2 = sortInd(sortInd_coin_cleaned+1);

[~, iInd_coin1, ~] = intersect(Ind_coin1, Ind_511);
[~, iInd_coin2, ~] = intersect(Ind_coin2, Ind_511);

iInd_coin = intersect(iInd_coin1,iInd_coin2);
Ind_coin1_511 = Ind_coin1(iInd_coin);
Ind_coin2_511 = Ind_coin2(iInd_coin);

Ind_coin = [Ind_coin1_511 Ind_coin2_511];


