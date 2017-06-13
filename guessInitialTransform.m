function [tform,gradError,xfit,zfit,xdat,zdat,xRed,zRed,stpR,stpG] = guessInitialTransform(redIm,greenIm)

[xRed,zRed] = FastPeakFind(redIm);
[xGreen,zGreen] = FastPeakFind(greenIm);

for m=1:length(xRed)
    lcs = [1:m-1,m+1:length(xRed)]; lcs = lcs(lcs>=1);
    r2 = sqrt((xRed(m)-xRed(lcs)).^2+(zRed(m)-zRed(lcs)).^2);
    stpR(m) = min(r2);
end
for m=1:length(xGreen)
    lcs = [1:m-1,m+1:length(xGreen)]; lcs = lcs(lcs>=1);
    r2 = sqrt((xGreen(m)-xGreen(lcs)).^2+(zGreen(m)-zGreen(lcs)).^2);
    stpG(m) = min(r2);
end

figure; plot(stpR); hold on; plot(stpG); hold off; drawnow;
    

r1 = sqrt((xGreen(1)-xRed).^2+(zGreen(1)-zRed).^2);
[mn,lc] = min(r1);
xred1 = xRed(lc); zred1 = zRed(lc);

r2 = sqrt((xGreen(end)-xRed).^2+(zGreen(end)-zRed).^2);
[mn,lc] = min(r2);
xred2 = xRed(lc); zred2 = zRed(lc);

aGuess = (xGreen(end)-xGreen(1))./(xred2-xred1);
eGuess = xGreen(1)-xred1.*aGuess;
fGuess = zGreen(1)-zred1.*aGuess;

figure;
startGuess = [aGuess,eGuess,fGuess];
model = @fun;
options = optimset('MaxFunEvals',500*length(startGuess),'MaxIter',500*length(startGuess));
est = fminsearch(model,startGuess,options);
[sse,xfit,zfit,xdat,zdat,gradError] = model(est);

figure; plot(xfit,zfit);

tform = affine2d([est(1),0,0;0,est(1),0;est(2),est(3),1]);

function [sse,xfit,zfit,xdat,zdat,errorSlope] = fun(params)
    a = params(1);
    e = params(2);
    f = params(3);
    
    xfit = xRed.*a+e;
    zfit = zRed.*a+f;
    
    for i=1:length(xfit)
        r = sqrt((xfit(i)-xGreen).^2+(zfit(i)-zGreen).^2);
        [mn,lc] = min(r);
        mins(i) = mn;
        xdat(i) = xGreen(lc);
        zdat(i) = zGreen(lc);
    end
    scatter(xGreen,zGreen); hold on; scatter(xfit,zfit); hold off; drawnow;
    sse = nansum(mins);
    
%     RMSEx = sqrt(sum((xfit-xdat).^2)/(length(xfit)-2));
%     RMSEz = sqrt(sum((zfit-zdat).^2)/(length(zfit)-2));
%     errorSlopeX = sqrt(sum((xfit-xdat).^2)/(length(xfit)-2))/sqrt(sum((xRed-mean(xRed)).^2))
%     errorSlopeZ = sqrt(sum((zfit-zdat).^2)/(length(zfit)-2))/sqrt(sum((zRed-mean(zRed)).^2))
%     errorSlope = sqrt(errorSlopeX.^2+errorSlopeZ.^2)
     
    ciX = confint(fit(xRed.',xdat.','poly1'));
    X = coeffvalues(fit(xRed.',xdat.','poly1'));
    ciZ = confint(fit(zRed.',zdat.','poly1'));
    Z = coeffvalues(fit(zRed.',zdat.','poly1'));
    erX = (abs(ciX(1,1)-X(1))+abs(ciX(2,1)-X(1)))/2;
    erZ = (abs(ciZ(1,1)-Z(1))+abs(ciZ(2,1)-Z(1)))/2;
    
    errorSlope = sqrt(erX.^2+erZ.^2);
    
   
    
end


end
    
    
    
    