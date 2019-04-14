function [icoefs,inums] = findpeakscwt(idata,scales,wavelet,idxStart,fplot)

    if nargin < 5
        % uncomment this to plot scalogram
        fplot = [];
    end

    if nargin < 4
        idxStart = 1;
    end
    bgSkip = 20;

%     ipeaks = zeros(size(idata));
    icoefs = zeros(size(idata)); % sum of coeffcients
    %inums = zeros(size(idata)); % sum of peak numbers on scalogram
   
    cwtData = cwt(idata,scales,wavelet);
    inums = zeros(length(scales),1);

    for r = 1:size(cwtData,1)
        cwtRow = cwtData(r,:);
        bgmax = 0;
%         take the first pixels till ysrf to be the background
        if idxStart <= bgSkip || idxStart - bgSkip > 3600
            bgmax = 0;
        else
            bgSig = cwtRow(1:idxStart-bgSkip);
            % maximum as threshold
            bgmax = max(bgSig);
        end
        
        % findpeaks along the cwt row
        try
        [pks,locs] = findpeaks(cwtRow);%,'minpeakheight',bgmax);
        locs(pks <= bgmax) = []; 
        catch
            bgmax
        end
        % if there is no peaks found in this cwt row, continue in next row
        locs(locs < idxStart) = [];
        if isempty(locs)
            continue;
        end
        
        % move the pixels to the location of maximum
        if ~strcmp(fplot,'plot')
            
            halfScale = floor(scales(r)/2);
            for j = 1:length(locs)
                if locs(j) + halfScale > 3600
                    continue;
                end
                sigSegment = idata(locs(j) - halfScale: locs(j) + halfScale);
                [~,ind] = max(sigSegment);
                if locs(j) - halfScale - 1  + ind < 3600
                    locs(j) = locs(j) - halfScale - 1  + ind;
                end
            end
        end

        icoefs(locs) = icoefs(locs) + cwtRow(locs)';
        %inums(locs) = inums(locs) + 1;
        inums(r) = length(locs);
        
        
        
% uncomment this to plot the local maxima on scalogram
        if strcmp(fplot,'plot')
            yyaxis left;
            bar(inums,'edgecolor','c','facecolor','c');
            ylabel('Nos. of CWT Maxima');
            set(gca,'YColor','k','Ydir','normal');
            
            yyaxis right;
            ylim([1 23]);
            if strcmp(wavelet,'mexh')
                colorline = 'b.';
            else
                colorline = 'b.';
            end
            plot(locs - 1,ones(size(locs)) * r,colorline);
            scaletitle = ['scales (1-' num2str(scales(end)) ')']; 
            ylabel(scaletitle);
            set(gca,'YColor','k');
        end

    end
end
