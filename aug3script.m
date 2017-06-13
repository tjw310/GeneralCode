function [X,Z] = aug3script(images,fixed,dyn,magCal)

[x,z]=FastPeakFind(images.projections(:,:,1).',1);

k=1;
for i=1:length(fixed)
    r = sqrt((fixed(i).centreX-x).^2+(fixed(i).centreY-z).^2);
    [mn,lc] = min(r);
    if mn<100
        xe(k) = x(lc);
        ze(k) = z(lc);
        xp(k) = fixed(i).centreX;
        zp(k) = fixed(i).centreY;
        y(k) = fixed(i).zDepth;
        k=k+1;
    end
end

figure; magCal.scatterAllPoints(dyn,'dynamic'); 
hold on; scatter(xe,ze,'gx'); scatter(xp,zp,'wx'); hold off;

e = findMagChanges(xp,xe,zp,ze,y,magCal);

x0 = magCal.opticalCentre(1);
z0 = magCal.opticalCentre(2);

trueXp = (xp-x0).*(1-e(1).*y)+x0;
trueXe = (xe-x0).*(1-e(2).*y)+x0;

X = (trueXp+trueXe)/2;

trueZp = (zp-z0).*(1-e(1).*y)+z0;
trueZe = (ze-z0).*(1-e(2).*y)+z0;

Z = (trueZp+trueZe)/2;

xp2 = (X-x0)./(1-e(1).*y)+x0;
zp2 = (Z-z0)./(1-e(1).*y)+z0;

xe2 = (X-x0)./(1-e(2).*y)+x0;
ze2 = (Z-z0)./(1-e(2).*y)+z0;

figure; magCal.scatterAllPoints(fixed,'fixed'); hold on; scatter(xe,ze,'gx');
scatter(xp2,zp2,'wx'); hold off;

figure; imagesc(images.projections(:,:,1).'); hold on; scatter(xe,ze,'ro');
scatter(xp,zp,'gx'); scatter(xe2,ze2,'wx'); hold off;

dr1 = sqrt((xp2-xp).^2+(zp2-zp).^2);
dr2 = sqrt((xe2-xe).^2+(ze2-ze).^2);

figure; plot(dr1); hold on; plot(dr2); hold off;


end


