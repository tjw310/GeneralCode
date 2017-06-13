function [ estimates,refinemodel,flag,fit ] = fitLinear( xdata,ydata )

xin = xdata(and(~isnan(xdata),~isnan(ydata)));
yin = ydata(and(~isnan(xdata),~isnan(ydata)));

start_guesses = [(yin(end)-yin(1))/(xin(end)-xin(1)),min(yin)];
model = @fun;
refinemodel = @refine;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);
[estimates,~,flag] = fminsearch(refinemodel,estimates,options);

% figure;
%drawnow;

[~,fit] = model(estimates);
%figure;
%scatter(xin,yin,[],1:length(xin),'filled'); hold on; plot(xin,fit,'r'); hold off; drawnow;

%figure;
%scatter(yin,xin,[],yin,'filled'); hold on; plot(fit,xin,'r'); hold off; drawnow;

    function [sse,fit] = fun(params)
        A = params(1);
        B = params(2);
        
        fit = A.*xin+B;
        
        sse = sum((fit-yin).^2);
    end

    function [sse,fit] = refine(params)
        A = params(1);
        B = params(2);
        
        fit = A.*xin+B;
        
        dif = fit-yin;
        err = std(dif);
        
        dif = dif(abs(dif-mean(dif))<err);
        
        sse = nansum(dif.^2);
    end


end
