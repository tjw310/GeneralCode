function [mtz0,motorMotion,AoRoffset] = findTmagandMotorMotionv2(obj)
xdata = obj.peaksXsub;

ymot = 1.3330*(obj.piezomotion-mean(obj.piezomotion));
xshift = circshift(ymot,[0,100]);

start_guesses = [obj.mtz0,20,0];
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));

figure;
estimates= fminsearch(model,start_guesses,options);
[~,motorMotion] = model(estimates);

mtz0 = estimates(1)
AoRoffset = estimates(2)
%plot(motorMotion); ylabel('Mean Motor Motion'); xlabel('Proj No'); title(sprintf('With AoRoffset = %.2f',AoRoffset));

    function [sse,trueDf] = fun(params)
            mt = params(1);
            shift = params(2);

            trueX = xshift./obj.pxSz.*mt+shift;
            xfit = (trueX-obj.opticCentre(1))./(1-obj.dmtdy*ymot)+obj.opticCentre(1);
            df = abs(xdata-xfit);
            df(df>200)=NaN;
            sse = nansum(df.^2);
            
            plot(xdata); hold on; plot(xfit); hold off; drawnow;

            dfScale = xdata-xfit;
            trueDf = dfScale.*(1-obj.dmtdy*ymot);
    end
end