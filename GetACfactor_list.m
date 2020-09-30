function [cilist,gooddetind] = GetACfactor_list(ci, detid_pair, numdet, R1, distrange)

[~, indang, ~, inddist, gooddetind] = list2sino(detid_pair, numdet, R1, distrange);

sz = [max(inddist), max(indang)];
if(any(sz-size(ci)))
    error('Dimension does not match!')
end

Ind = sub2ind(sz,inddist,indang);

cilist = NaN(size(detid_pair,1),1);
cilist(gooddetind) = ci(Ind);
