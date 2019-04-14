% function peakim = peakimcwt(im,scales,wavelet,ysrf,ybtm,app,bgSkip)
function [peakim,nums] = peakimcwt(im,scales,wavelet,ysrf)
%--------------------------------------------------------------------------
% Generate the feat image by using average of CWT coefficient
%
% Parameters
% ----------
% im : 2-D image
%   SHARAD echogram
% wavelet : string
%   Name of mother wavelet, 'mexh' or 'mexh'
% scales : 1 * n array
%   Withs of the mother wavelet, [1,4,8,16]
% bgSkip : int
%   Pixels below bed reflectors to start calculating background noise, 50
% 
% See also
% --------
% read_data.m
% preprocessing.m
%
% Returns
% -------
% peakim : 2-D image
%   Feature image as same dimension as 'im'
%   Pixel values of peaks are wavelet coefficient
%   Non-peaks are indicated by zero.
%
% Examples
% --------
% [geoinfo,echogram] = read_data(inputfilename);
% [im_nl,im_db, ysrf,ybtm] = preprocessing(geoinfo,echogram);
% peakim = find_peaks_cwt(im,ysrf,ybtm);
%--------------------------------------------------------------------------

    %% Define parameters

    if nargin < 3
        wavelet = 'mexh';
    end
    
    if nargin < 2
        scales = 1:16;
    end

    peakim = zeros(size(im));
    nums = zeros(length(scales), size(im,2));

    %% Start processing
    for i = 1: size(im,2)
        [peakim(:,i), nums(:,i)] = findpeakscwt(im(:,i), scales,wavelet,ysrf(i));
    end

end

% % % function peaks = SHARADCWTfeatim(im,ysrf,SCALE_RANGE, WAVELET_TYPE,BG_SKIP)
% % % [~, num_bins] = size(im);
% % % peaks = zeros(size(im));
% % % for i = 1:num_bins
% % %     iprofile = im(:,i);
% % %     CWTcoeffs = cwt(iprofile,SCALE_RANGE,WAVELET_TYPE); % colormap jet;
% % %     mCWTcoeffs = mean(CWTcoeffs);
% % %     if ysrf(i) - BG_SKIP > 0 
% % %         bgnoise_top = mCWTcoeffs(1:ysrf(i) -BG_SKIP);
% % %         thresh_top = max(bgnoise_top);
% % %     else
% % %         thresh_top = -9999;
% % %     end
% % %     bgnoise_btm = mCWTcoeffs(3600-BG_SKIP:3600);
% % %     thresh = max(max(bgnoise_btm),thresh_top);
% % %         
% % %     [pks,locs] = findpeaks(mCWTcoeffs,'minpeakheight',0.8*thresh);
% % %     peaks(locs,i) = pks;
% % %     peaks(1:ysrf(i),i) = 0;
% % %     peaks(3600-BG_SKIP:3600,i) = 0;
% % % end
% % % 
% % % end