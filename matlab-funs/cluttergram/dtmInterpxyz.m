%% General utility that can interpolate DTM values according to x, y coordinates
% Si-Ting Xiong
% Compiled on 11-03-2017


function vdtm = dtmInterpxyz(dtm,dtminfo,x,y,method)
ulx = dtminfo(1);
uly = dtminfo(2);
resx = dtminfo(3);
resy = dtminfo(4);

% [nx, ny] = size(dtm);
% vdtm = zeros(length(x),1);
% for i = 1:length(x)
%     xi = floor((x(i) - ulx) / resx);
%     yi = floor((y(i) - uly) / resy);
%     
%     if xi > 0 && xi <= nx && yi > 0 && yi <= ny
%         vdtm(i) = dtm(xi,yi) + 3396000;
%     else
%         vdtm(i) = nan;
%     end
% end
[nx,ny] = size(dtm);
minX = min(x); maxX = max(x);
minY = min(y); maxY = max(y);

istart = floor((minX - ulx) / resx) - 5;
iend = ceil((maxX - ulx) / resx) + 5;
jend = floor((minY - uly) / resy) + 5;
jstart = ceil((maxY - uly) / resy) - 5;

XX = ulx:resx:ulx + (resx * (ny-1));
YY = uly:resy:uly + (resy * (nx-1));

%%% ADD WARNING HERE IF THE EXTEND out DTM
if istart > 0 && iend <= ny && jstart > 0 && jend <= nx
    dtm = dtm(jstart:jend,istart:iend);
    xx = XX(istart:iend);
    yy = YY(jstart:jend);
else
    xx = XX;
    yy = YY;
end

[X,Y] = meshgrid(xx,yy);
% surf(X,Y,subdtm);

vdtm = interp2(X,Y,dtm,x,y,'nearest') + 3396000;

end
