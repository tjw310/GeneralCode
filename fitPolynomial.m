function [ estimates,refinemodel,flag ] = fitPolynomial( xdata,ydata )

start_guesses = [1,1,1,1,1,(ydata(end)-ydata(1))/(xdata(end)-xdata(1)),min(ydata)];
model = @fun;
refinemodel = @refine;
options = optimset('MaxFunEvals',2000*length(start_guesses),'MaxIter',2000*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);
[estimates,~,flag] = fminsearch(refinemodel,estimates,options);

% figure;
%drawnow;

    function [sse,fit] = fun(params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        E = params(5);
        F = params(6);
        G = params(7);
        
        fit = A.*xdata.^6 + B*xdata.^5 + C.*xdata.^4+D.*xdata.^3+E.*xdata.^2+F.*xdata+G;
        
        sse = sum((fit-ydata).^2);
    end

    function [sse,fit] = refine(params)
        A = params(1);
        B = params(2);
        C = params(3);
        D = params(4);
        E = params(5);
        F = params(6);
        G = params(7);
        
        fit = A.*xdata.^6 + B*xdata.^5 + C.*xdata.^4+D.*xdata.^3+E.*xdata.^2+F.*xdata+G;
        
        dif = fit-ydata;
        err = 2*std(dif);
        
        %disregard top 10% of points.
        d = abs(dif);
        d = sort(d,'descend');
        d = d(round(length(d)/10):end);
        
        %dif = dif(abs(dif-mean(dif))<err);
        
        sse = sum(d.^2);
    end


end
