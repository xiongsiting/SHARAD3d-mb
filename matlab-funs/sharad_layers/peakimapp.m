function peakim = peakimapp(im,ysrf)
%-----------------------------------------------------
% Generate the feat image by using Raquel's APP method
%
% Parameters
% ----------
% im : 2-D image
%   MCoRDS echogram
% ysrf : 1 * n array
%   Row numbers of surface echoes
% ybtm : 1 * n array
%   Row numbers of bottom echoes
% lwins : 1 * n array
%   A set sizes for long windows
% swins : 1 * n array
%   A set sizes for short windows
%   same size length with lwins.
% weights : 1 * n array
%   A set sizes for weights
%   same size with lwins.
%
% Returns
% -------
% peakim : 2-D image
%   Feature image as same dimension as 'im'
%   Pixel values of peaks are STA/LTA ratio
%   Non-peaks are indicated by zero.
%
% Examples
% --------
% [geoinfo,echogram] = read_data(inputfilename);
% [im_nl,im_db, ysrf,ybtm] = preprocessing(geoinfo,echogram);
% peakim_eda = find_peaks_eda(im_nl,ysrf,ybtm);
%-----------------------------------------------------
    if nargin < 3
        ybtm = ones(1,size(im,2)) * size(im,1);
    end
    
    if nargin < 2
        ysrf = ones(1,size(im,2));
    end
    
    [nrows,ncols] = size(im);
    peakim = zeros(nrows,ncols);
    
    for i = 1:ncols
        % initial estimation using the first set of parapmeters
        peakim(:,i) = findpeaksapp(im(:,i),ysrf(i));
    end
%     peakim = scale2minmax(peakim);
%     peakim_nan = peakim;
%     peakim_nan(peakim == 0) = nan;
%     minval = min(peakim_nan(:));
%     maxval = max(peakim_nan(:));
%     peakim = (peakim_nan - minval)/(maxval-minval);
%     peakim(isnan(peakim)) = 0;
end
