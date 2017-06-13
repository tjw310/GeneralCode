function estimates = findTmagv2(obj,xdata)

ymot = obj.piezomotion;

start_guesses = [obj.mtz0];
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);

function sse = fun(params)
        mt = params(1);
        xfit = circshift(ymot,[100,0])./obj.pxSz.*mt;
        
        
        xfit = (xfit-obj.opticCentre(1))./(1-obj.dmtdy*ymot)+obj.opticCentre(1);
        
        sse = sum((xfit-xdata).^2);
        
        scatter(1:400,xdata); hold on; plot(xfit); hold off; pause(1); drawnow;
        
end

end