function [M,dose_data,masks]=BuildDoseMatrix_CTsim(h5file, maskfile, fmapsfile)

verbose = 2;
PTVind = 1;

[ masks ] = open_masks( maskfile, 'xyz', verbose );
PTV = masks{PTVind}.mask;
sz = size(PTV);

for i = 1:length(masks)
    masks{i}.mask = permute(masks{i}.mask,[2,1,3]);
    masks{i}.name = masks{i}.name;
end

M = read_sparse_mcdose(h5file);

dose_data = read_fmaps(fmapsfile, verbose, 'all');

[r,c,v] = find(M);
[Nrows, Ncols] = size(M);
dicomsize = sz;
[I,J,K] = ind2sub(dicomsize,r);
dicomsizenew = dicomsize;
dicomsizenew(2) = dicomsize(1);
dicomsizenew(1) = dicomsize(2);
r = sub2ind(dicomsizenew,J,I,K);

M = sparse(r, c, v, Nrows, Ncols); %OUTVAR