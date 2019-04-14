function [ipeaks,ikurto] = findpeaksapp(idata,idxStart)
    
    if nargin < 5
        idxEnd = length(idata);
    end
    
    if nargin < 2
        idxStart = 1;
    end
    
    weights = [0.3,0.25,0.2,0.15,0.05,0.05];
    lwins = [5,7,9,11,13,13];
    kSigma = 3;
    swins = [3,5,5,7,7,9];
        
    n = length(idata);
    ipeaks = zeros(size(idata));
    ikurto = zeros(size(idata));
    
    % check for local maxima
    [~,relmax] = findpeaks(idata);
    
    relmax(relmax <= idxStart) = [];
    relmax(relmax >= idxEnd) = [];
    
    for k = 1:length(swins)
        % half window sizes
        try
        hlwin = int16(floor(lwins(k)/2));
        hswin = int16(floor(swins(k)/2));
        catch
            hlwin
        end
        
        for j = 1:length(relmax)
            i = relmax(j); 
            % deal with head and end part
            if i <= hlwin
                sterm = idata(hlwin + 1 - hswin: hlwin + 1 + hswin);
                lterm = idata(1:lwins(k));
            else
                if i > n - hlwin
                    sterm = idata(n - hlwin - hswin: n - hlwin + hswin);
                    lterm = idata(n - lwins(k):end);
                else
                    sterm = idata(i-hswin:i+hswin);
                    lterm = idata(i-hlwin:i+hlwin);
                end
            end
            % calculate the STA/LTA ratio, not exactly average ???
            ipeaks(i) = ipeaks(i) + weights(k) * (sum(sterm)/sum(lterm));
            ikurto(i) = ikurto(i) + weights(k) * kurtosis(lterm);
        end
    end
    if kSigma > 0
        idx = find(ipeaks * kSigma < ikurto);
        ipeaks(idx) = 0;
        ikurto(idx) = 0;
    end
    
    bgSig = ipeaks(1:idxStart);
    % maximum as threshold
    bgmax = max(bgSig);
    ipeaks(ipeaks < bgmax) = 0;
    

end
