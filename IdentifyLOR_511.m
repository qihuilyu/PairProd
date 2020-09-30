function Ind_coin = IdentifyLOR_511(energy, sortedtime, sortInd, CoincidenceTime)

Ind_511 = find(abs(energy-0.511)<0.0001);

sortInd_coin = find(diff(sortedtime)<CoincidenceTime);
Ind_coin1 = sortInd(sortInd_coin);
Ind_coin2 = sortInd(sortInd_coin+1);

[~, iInd_coin1, ~] = intersect(Ind_coin1, Ind_511);
[~, iInd_coin2, ~] = intersect(Ind_coin2, Ind_511);

iInd_coin = intersect(iInd_coin1,iInd_coin2);
Ind_coin1_511 = Ind_coin1(iInd_coin);
Ind_coin2_511 = Ind_coin2(iInd_coin);

Ind_coin = [Ind_coin1_511 Ind_coin2_511];

