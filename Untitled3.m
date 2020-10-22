
[C,ia,ic] = unique(detectorIds(1:1000000));

counts = hist(detectorIds(1:1000000),1440);


B = groupcounts(detectorIds(ceil(rand(1000000,1)*numel(detectorIds))));
figure; plot(B)


B = groupcounts(detid_pair(ceil(rand(1000000,1)*numel(detid_pair(:)))));
figure; plot(B)


B = groupcounts(inddist(ceil(rand(1000000,1)*numel(inddist(:)))));
figure; plot(B,'o')


CorrectedTime = globalTimes + AlleventID*10;


B = groupcounts(detectorIds(ceil(rand(1000000,1)*numel(detectorIds))));
figure; plot(B)


B = groupcounts(detid_pair(ceil(rand(1000000,1)*numel(detid_pair(:)))));
figure; plot(B)


