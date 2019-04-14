function [Ps,surfind,surfhei] = cluttergram(dtmpatch,X,Y,geoinfo, mtype)
% Updated on 13-10-2017 adapted from previous version for package of SPLD-PL

% We post radargram columns at 128 points per degree spacing (about 460 m along track), 
% with the center (1800th) range cell corresponding to the free- space, round-trip time delay 
% of the MOLA-defined Mars areoid. 

% Simulate the cluttergram
% See from Adamo Ferro IEEE TGRS, VOL. 51 NO. 5 2013
C = 299792458;sigma = 2.6;tspace = 0.0375e-6;
ncols = size(geoinfo,1);
x = geoinfo(:,3);
y = geoinfo(:,4);
Hsat = geoinfo(:,6);
Ha = geoinfo(:,5);
time = 2*(Hsat - Ha)/C;

Ps = zeros(3600,ncols);
surfind = zeros(1,ncols);
% If pixels locate in Ai
for i = 1:ncols
    for j = 2:size(dtmpatch,1)-1
        X0 = x(i); Y0 = y(i); h0 = Hsat(i);
        Xp = X(j,i);Yp = Y(j,i); hp = dtmpatch(j,i);
            if isnan(hp)
            continue;
            end
        Rxy = sqrt((Xp - X0).^2 + (Yp - Y0).^2 + (hp-h0).^2);
        t = 2 * Rxy / C;
        tdelay = 1800 + round((t - time(i)) / tspace);
        if j == ceil(size(dtmpatch,1)/2)
            surfind(i) = tdelay;
        end
        switch mtype
        case '1'
            try
                Ps(tdelay,i) = Ps(tdelay,i) + 1/power(Rxy,4);
            catch
            i = i;
            end
        case '2'
            if j < ceil(size(dtmpatch,1)/2)
                %Using local incidence angle and relative dielectric constant
                pxspace = sqrt((X(j-1,i)-X(j,i)).^2 + (Y(j-1,i)-Y(j,i)).^2);
                % Local incidence angle
                alpha = atan(abs(dtmpatch(j-1,i)-dtmpatch(j,i))/pxspace);
                theta = acos((h0-hp)/Rxy) - alpha;
            else
            %Using local incidence angle and relative dielectric constant
            pxspace = sqrt((X(j+1,i)-X(j,i)).^2 + (Y(j+1,i)-Y(j,i)).^2);
            % Local incidence angle
            alpha = atan(abs(dtmpatch(j+1,i)-dtmpatch(j,i))/pxspace);
            theta = acos((h0-hp)/Rxy) - alpha;
            end
            % rho is the relative Fresnel coefficient of the surface
            rho = ((cos(theta)-sqrt(sigma - sin(theta)^2))/(cos(theta) + sqrt(sigma - sin(theta)^2))).^2;
            Ps(tdelay,i) = Ps(tdelay,i) + rho * cos(theta)/power(Rxy,4);
        case '3'
            %% The equation is (4) from M.G. Spagnuolo et al. Planetary and Space Science 59(2011) 1222-1230
            K = sqrt(0.01);
            pxspace = sqrt((X(j+1,i)-X(j,i)).^2 + (Y(j+1,i)-Y(j,i)).^2);
            if j < ceil(size(dtmpatch,1)/2)
                %Using local incidence angle and relative dielectric constant
                pxspace = sqrt((X(j-1,i)-X(j,i)).^2 + (Y(j-1,i)-Y(j,i)).^2);
                % Local incidence angle
                alpha = atan(abs(dtmpatch(j-1,i)-dtmpatch(j,i))/pxspace);
                theta = acos((h0-hp)/Rxy) - alpha;
            else
            %Using local incidence angle and relative dielectric constant
            pxspace = sqrt((X(j+1,i)-X(j,i)).^2 + (Y(j+1,i)-Y(j,i)).^2);
            % Local incidence angle
            alpha = atan(abs(dtmpatch(j+1,i)-dtmpatch(j,i))/pxspace);
            theta = acos((h0-hp)/Rxy) - alpha;
            end
            rho = ((cos(theta)-sqrt(sigma - sin(theta)^2))/(cos(theta) + sqrt(sigma - sin(theta)^2))).^2;
            rho = rho * K / 2 * power(power(cos(theta),4) + K * sin(theta)^2,-3/2);
            Ps(tdelay,i) = Ps(tdelay,i) + rho * cos(theta)/power(Rxy,4);
        end
    end
end
surfhei = Ha' + (1800 - surfind)*tspace*C/2;

end
