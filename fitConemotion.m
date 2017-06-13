function [ estimates,flag ] = fitConemotion( xdata,zdata,t,coneReconClass )

opticCentre = coneReconClass.opticCentre;

xdata = xdata - 1024;
zdata = zdata - 1024;

r = (max(xdata)-min(xdata))/2;
yGuess = sqrt(r^2-(xdata(1))^2);

start_guesses = [coneReconClass.AoRangle,10*coneReconClass.dmdP,xdata(1),0,zdata(1)]
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);

% figure;
%drawnow;

    function [sse,fit] = fun(params)
        x = params(3);
        y = params(4);
        z = params(5);
        %x = xdata(1);
        %y = yGuess;
        %z = zdata(1);
        %AoRcentre = coneReconClass.AoRcentre;
        AoRangle = params(1);
        %AoRangle = coneReconClass.AoRangle;
        dmdy = params(2);
        

        for i=1:length(t)
            theta = t(i)/180*pi;
            w = (1+tan(AoRangle/180*pi)^2)^-0.5;
            u = w*tan(AoRangle/180*pi);
            v=0;
        
        x1 = u.*(u.*x+v.*y+w.*z)*(1-cos(theta))+x*cos(theta)+(-w.*y+v.*z).*sin(theta);
        y1 = v.*(u.*x+v.*y+w.*z)*(1-cos(theta))+y*cos(theta)+(w.*x-u.*z).*sin(theta);
        z1 = w.*(u.*x+v.*y+w.*z)*(1-cos(theta))+z*cos(theta)+(-v.*x+u.*y).*sin(theta);
        

        x2(i) = (x1-opticCentre(1))./(1-dmdy.*y1)+opticCentre(1);
        z2(i) = (z1-opticCentre(2))./(1-dmdy.*y1)+opticCentre(2);
        
        end
        
        scatter(xdata,zdata); hold on; scatter(x2,z2); hold off; drawnow;
        xcopy = xdata; zcopy = zdata;
        sse = 0;
%         while ~isempty(x2)
%             dif = sqrt((xcopy-x2(1)).^2+(zcopy-z2(1)).^2);
%             [mn,lc]=min(dif);
%             xcopy = xcopy((1:length(xcopy))~=lc);
%             x2 = x2(2:end);
%             zcopy = zcopy((1:length(zcopy))~=lc);
%             z2 = z2(2:end);
%             sse = sse+mn;
%         end

       sse = sum(sqrt((xdata-x2).^2+(zdata-z2).^6));
        
        
    end

end
