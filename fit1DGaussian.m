function [ estimates,model,flag ] = fit1DGaussian( xdata,ydata )

[~,lc] = max(ydata);

sigGuess = 0.02;

start_guesses = [max(ydata) sigGuess xdata(lc) min(ydata) ]
    %mean(min(ydata(:)))+std(ydata(:))];

model = @fun;
options = optimset('MaxFunEvals',1e4*length(start_guesses),'MaxIter',1e4*length(start_guesses),'Tolx',1e-5,'tolfun',1e-5);
[estimates,~,flag] = fminsearch(model,start_guesses,options);

if flag~=1
    estimates(1:4) = NaN;
end
% figure;
%drawnow;
[~,fit]=model(estimates);


        plot(xdata,ydata);
        hold on
        plot(xdata,fit); hold off;
        title(num2str(flag));
        drawnow;

    function [sse,fit] = fun(params)
        A = params(1);
        sigx = params(2);
        x0 = params(3);
        B = params(4);
        
        %xdata = linspace(-1,1,size(ydata,2)).*size(ydata,2)/2;
        
        fit = A.*exp(-(xdata-x0).^2./(2*sigx^2))+B;
% %         
%         plot(xdata,ydata);
%         hold on
%         plot(xdata,fit); hold off;
%         drawnow;
% %         
        sse = sum(sum((fit-ydata).^2));
    end


end
