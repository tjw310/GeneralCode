function [mtz0,motorMotion,AoRoffset] = findTmagandMotorMotion(obj,xdata)

ymot = 1.3330*(obj.piezomotion-mean(obj.piezomotion));

start_guesses = [obj.mtz0,0];
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses),'TolX',1e-8,'TolFun',1e-8);

figure;
disp('run')
estimates= fminsearch(model,start_guesses,options);
[~,motorMotion] = model(estimates);
disp('run');
mtz0 = estimates(1);
AoRoffset = estimates(2);
%plot(motorMotion); ylabel('Mean Motor Motion'); xlabel('Proj No'); title(sprintf('With AoRoffset = %.2f',AoRoffset));

    function [sse,trueDf] = fun(params)
            mt = params(1);
            AoRoffset = params(2);

            trueX = circshift(ymot,[0,100])./obj.pxSz.*mt+AoRoffset;
            xfit = (trueX-obj.opticCentre(1))./(1-obj.dmtdy*ymot)+obj.opticCentre(1);
            sse = sum((xdata-xfit).^2);
            
            plot(xdata); hold on; plot(xfit); hold off; drawnow; pause(1);

            dfScale = xdata-xfit;
            trueDf = dfScale.*(1-obj.dmtdy*ymot);
    end
end