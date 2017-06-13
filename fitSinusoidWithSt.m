function [ estimates,fit,flag ] = fitSinusoidWithSt(ydata,start_guesses)


if size(ydata,1)~=1
    ydata = ydata.';
end

model = @fun;
options = optimset('MaxFunEvals',5000*length(start_guesses),'MaxIter',5000*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);

model = @refine;
options = optimset('MaxFunEvals',5000*length(start_guesses),'MaxIter',5000*length(start_guesses));
[estimates,~,flag] = fminsearch(model,estimates,options);

% figure;
%drawnow;

[~,fit] = model(estimates);
  %figure;
  %plot(ydata); hold on; plot(fit,'r'); hold off; drawnow;

    function [sse,fit] = fun(params)
        A = params(1);
        B = params(2);
        C = params(3);
        
        theta = linspace(0,2*pi-2*pi/length(ydata),length(ydata));
        
        fit = A.*sin(theta+B)+C;
        
        sse = nansum((fit-ydata).^2);
    end

    function [sse,fit] = refine(params)
        A = params(1);
        B = params(2);
        C = params(3);
        
        theta = linspace(0,2*pi-2*pi/length(ydata),length(ydata));
        
        fit = A.*sin(theta+B)+C;
        
        df = abs(fit-ydata);
        df(df>100)=NaN;
        
        sse = nansum(df.^2);
    end

end
