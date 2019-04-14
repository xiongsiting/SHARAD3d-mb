function peakim = removeclutter(Ps,surfind,peakim,surfecho)

%surfclut = pickFirstReturn(Ps, 'mouginot');
surfclut = size(1,size(Ps,2));
for i = 1:size(Ps,2)
    segment = Ps(surfind(i) - 50:surfind(i) + 50,i);
    [~, maxind] = max(segment);
    surfclut(i) = surfind(i) - 50 - 1 + maxind;
end

clutim = zeros(size(peakim));
for i = 1:size(clutim,2)
    bgnoise = max(Ps(surfclut(i)-10:surfclut(i) + 10,i));
    [~,locs] = findpeaks(Ps(:,i),'minpeakheight',bgnoise * 0.5);

    clutim(locs,i)= 1;
end
%% subsurface
% % delth = round(median(surfecho -surfclut));
% % if abs(delth) <= 15 && delth ~= 0
% %     for i = 1:size(peakim,2)
% %         pkind = find(peakim(:,i) > 0);
% %         pkind(pkind > 3600 + delth) = [];
% %         if isempty(pkind)
% %             continue;
% %         end
% %         peakim(pkind - delth,i) = peakim(pkind,i);
% %         peakim(pkind,i) = 0;
% %     end
% % end

for i = 1:size(peakim,2)
    delth = surfecho(i) - surfclut(i);
    pkind = find(peakim(:,i) > 0);
    pkind(pkind > 3600 + delth) = [];
    if isempty(pkind) || delth == 0
        continue;
    else
        if delth > 15
            peakim(pkind,i) = 0;
        end
    end
    peakim(pkind - delth,i) = peakim(pkind,i);
    peakim(pkind,i) = 0;
end

for i = 1:size(peakim,2)
%     peakim(1:surfind(i)-1,i) = 0;

    try
        clind = find(clutim(:,i) > 0);
        pkind = find(peakim(:,i) > 0);
    catch
        i
    end
    
    if isempty(pkind) || isempty(clind)
        % comment following lines for removing surface layer
%         peakim(surfind(i),i) = 1;
        continue;
    end

    lind = max(1,pkind(1) - 1);
    uind = min(3600, pkind(1) + 1);
    peakim(lind:uind,i) = 0;

    for k = 1:length(clind)
        lind = max(1, clind(k) - 5);
        uind = min(3600, clind(k) + 5);
        peakim(lind:uind,i) = 0;
    end
    
    %peakim(surfind(i),i) = 1;
    peakim(surfind(i),i) = 0; % try to remove surface layer
end

end