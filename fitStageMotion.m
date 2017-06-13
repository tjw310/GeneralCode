function [ estimates,sinu,flag ] = fitStageMotion(ydata,stepsize)

start_guesses = fitSinusoid(ydata);
model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);

[~,fit] = model(estimates);
%figure;
%plot(ydata); hold on; plot(fit,'r'); hold off; drawnow;

    function [sse,fit] = fun(params)
        A = params(1);
        B = params(2);
        C = params(3);
        
        theta = linspace(0,2*pi-2*pi/length(ydata),length(ydata));
        
        sinu = A.*sin(theta+B)+C;
        fit = round(sinu/stepsize)*stepsize;
        sse = sum((fit-ydata).^2);        
        
        %plot(ydata); hold on; plot(sinu); plot(fit); hold off; drawnow;
        
    end

end
