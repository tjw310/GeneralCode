function estimates = fitAugust31(xOut,zOut,progX)  

x = xOut-repmat(mean(xOut,1),size(xOut,1),1);
z = zOut-repmat(mean(zOut,1),size(xOut,1),1);

r = sqrt(x.^2+z.^2);
theta = atan2(z,x);
r(theta>=pi/2) = -r(theta>=pi/2);

dr = r(2:end,:)-r(1:end-1,:);
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

    fit = repmat((A.*dp).',1,size(x,2));

    subplot(2,2,1); 
    for i=1:size(dr,2)
    scatter(1:size(dr,1),dr(:,i),'b'); hold on; plot(fit(:,i),'r'); 
    end; 
    hold off;
    subplot(2,2,2)
    plot(fit-dr); drawnow;
    
    x2 = A.*(progX-mean(progX));
    subplot(2,2,3)
    for i=1:size(dr,2)
    scatter(1:size(r,1),r(:,1),'b'); hold on;
    end; plot(x2,'r'); hold off;
    subplot(2,2,4)
    plot(repmat(x2.',1,size(dr,2))-r); drawnow;

    sse = sum(sum((fit-dr).^2));
end

end