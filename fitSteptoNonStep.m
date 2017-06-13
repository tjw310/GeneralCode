function [mag,offset,dfStep,motorMotion] = fitSteptoNonStep(obj1,obj2)
% object 2 has the stage motion
x1 = obj1.peaksXsub;
x2 = obj2.peaksXsub;

model = @fun;
model2=@refine;
st_guess = [obj2.mtz0,10];
estimates = fminsearch(model,st_guess);
estimates = fminsearch(model2,estimates);
[~,dfStep] = model2(estimates);
mag = estimates(1);
offset = estimates(2);

[e,f] = fitSinusoid(x1);
motorMotion = x1-f;
motorMotion(isnan(dfStep))=NaN;

[e2,f2] = fitSinusoid(mag.*obj2.stagemotion/obj2.pxSz+x2);
motorMotion2 = mag.*obj2.stagemotion/obj2.pxSz+x2-f2;
motorMotion2(isnan(dfStep))=NaN;

dfStep = fitLinearToNan(dfStep);
motorMotion = fitLinearToNan(motorMotion);

figure; plot(dfStep); title('Difference In Stepped to Non-Stepped'); xlabel('Proj No'); ylabel('Pixels');
figure; plot(motorMotion); title('Motor Motion'); xlabel('Proj No'); ylabel('Pixels');
hold on; plot(motorMotion2); hold off;
function [sse,df] = fun(params)
    a = params(1);
    b = params(2);
    
    fit = x2+obj2.stagemotion/obj2.pxSz*a+b;
    
    %plot(x1);
   % hold on; plot(fit); hold off; drawnow; pause(0.5);
    
    df = fit-x1;
    df(abs(df)>100) = NaN;
    sse = nansum(df.^2);
    
end

function [sse,df] = refine(params)
    a = params(1);
    b = params(2);
    
    fit = x2+obj2.stagemotion/obj2.pxSz*a+b;
    
    plot(x1);
    hold on; plot(fit); hold off; drawnow;
    
    df = fit-x1;
    
    df(abs(df)>50)=NaN;
    sse = nansum(df.^2);
    
end

end