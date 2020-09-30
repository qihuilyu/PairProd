function [out, sino_filt] = fbp2_QL(sino, geom, varargin)

if geom.sg.nb ~= size(sino, 1) || geom.sg.na ~= size(sino, 2)
    error 'bad sino size'
end

opt.window = '';
opt = vararg_pair(opt, varargin);

sino_filt = fbp2_sino_filter('flat', sino, ...
    'ds', geom.sg.dr, 'window', opt.window);

out = fbp2_back_QL(geom.arg_back2{:}, single(sino_filt));


