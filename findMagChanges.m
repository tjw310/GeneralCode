function [ estimates,flag,trueX,trueZ,xGuess,zGuess] = findMagChanges( xp,xe,zp,ze,y,magCal )

start_guesses = [-0.35,0.1];
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses),'tolx',1e-8,'tolfun',1e-8);
[estimates,~,flag] = fminsearch(model,start_guesses,options);

[~,trueX,trueZ,xGuess,zGuess,medDist]= model(estimates);
figure;
scatter(xe,ze,'b'); hold on; scatter(xGuess,zGuess,'rx'); hold off; drawnow; title(num2str(medDist)); drawnow;

figure;
scatter(xe,ze,[],y,'filled'); hold on; scatter(magCal.opticalCentre(1),magCal.opticalCentre(2),'k'); 
scatter(trueX,trueZ,'rx'); hold off;

    function [sse,x,z,xeg,zeg,sse2] = fun(params)
        dmpdy = params(2);
        %dmpdy = 0.1978;
        dmedy = params(1);
        %const = params(3);
        xo = magCal.opticalCentre(1);
        zo = magCal.opticalCentre(2);
        
       % xo = params(2);
        %zo = params(3);
        
        x = (xp-xo).*(1-dmpdy.*y)+xo;
        z = (zp-zo).*(1-dmpdy.*y)+zo;
        
        xeg = (x-xo)./(1-dmedy.*y)+xo;
        zeg = (z-zo)./(1-dmedy.*y)+zo;
        
        %y1 = (xp-xe)./(dmpdy.*(xp-x0)+dmedy.*(x0-xe));
        
        %y2 = (zp-ze)./(dmpdy.*(zp-z0)+dmedy.*(z0-ze));
        
        %plot(1:length(y),y1-y); hold on; plot(1:length(y),y2-y,'r'); hold off; drawnow;
        
        %sse = sum(sqrt((y1-y).^2+(y2-y).^2));
        
        sse = nansum(sqrt((xe-xeg).^2+(ze-zeg).^2));
        sse2 = nanmedian(sqrt((xe-xeg).^2+(ze-zeg).^2));
        
       % scatter(xe,ze,'b'); hold on; scatter(xeg,zeg,'rx'); hold off; drawnow;
        
    end

end
