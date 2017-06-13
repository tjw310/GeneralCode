function estimates = fitStageMotionToProgValues(xOut,zOut,progX)  

x = xOut-mean(xOut);
z = zOut-mean(zOut);

r = sqrt(x.^2+z.^2);
theta = atan2(z,x);
r(theta>=pi/2) = -r(theta>=pi/2);

dr = r(2:end)-r(1:end-1);
dp = progX(2:end)-progX(1:end-1);

st = [-.05];
model = @fitM;
figure;
options = optimset('tolx',1e-8,'tolfun',1e-8);
estimates = fminsearch(model,st);

sprintf('mag is %f',abs(6.5e-3*estimates(1)))

%[~,x2] = model(estimates);
%figure; plot(x2(1:100)-r(1:100)); hold on; plot(-1*(x2(101:end)-r(101:end))); hold off;


function [sse,x2] = fitM(params)
    A = params(1);

    fit = A.*dp;

    subplot(2,2,1)
    scatter(1:length(dr),dr); hold on; plot(fit); hold off;
    subplot(2,2,2)
    plot(fit-dr); drawnow;
    
    x2 = A.*(progX-mean(progX));
    subplot(2,2,3)
    scatter(1:length(r),r); hold on; plot(x2); hold off;
    subplot(2,2,4)
    plot(x2-r); drawnow;

    sse = sum((fit-dr).^2);
end

end