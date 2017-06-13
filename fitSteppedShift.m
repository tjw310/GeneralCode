function [mtz0,motorMotion,AoRoffset] = fitSteppedShift(obj,xdata)

ymot = 1.3330*(obj.piezomotion-mean(obj.piezomotion));

start_guesses = [obj.mtz0,2*mean(obj.stagemotion)/obj.pxSz*obj.mtz0];
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses),'TolX',1e-16,'TolFun',1e-16);

figure;
estimates = fminsearch(model,start_guesses,options);
[~,motorMotion] = model(estimates);

mtz0 = estimates(1);
AoRoffset = estimates(2);
figure;
plot(motorMotion); ylabel('Mean Motor Motion'); xlabel('Proj No'); title(sprintf('With AoRoffset = %.2f',AoRoffset));

    function [sse,trueDf] = fun(params)
            mt = params(1);
            offset = params(2);

            xproj = (circshift(ymot,[0,obj.nProj/4])-obj.stagemotion)/obj.pxSz.*mt+offset;      
            xfit = (xproj-obj.opticCentre(1))./(1-obj.dmtdy*ymot)+obj.opticCentre(1);
            sse = sum((xdata-xfit).^2);
            
            plot(xdata); hold on; plot(xfit); hold off; drawnow;

            dfScale = xdata-xfit;
            trueDf = dfScale.*(1-obj.dmtdy*ymot);
    end

end