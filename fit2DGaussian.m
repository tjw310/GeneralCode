function [ estimates,model,flag ] = fit2DGaussian( xydata )

[row,col] = find(xydata==max(xydata(:)));

row = round(mean(row));
col = round(mean(col));

sigGuess = size(xydata,2)/5;

start_guesses = [mean(max(xydata(:))) 0 sigGuess sigGuess col-size(xydata,2)/2 row-size(xydata,1)/2 mean(min(xydata(:)))+std(xydata(:))];

model = @fun;
options = optimset('MaxFunEvals',500*length(start_guesses),'MaxIter',500*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);

% figure;
%drawnow;

    function [sse,fit] = fun(params)
        A = params(1);
        theta = params(2);
        sigx = params(3);
        sigy = params(4);
        x0 = params(5);
        y0 = params(6);
        B = params(7);
        
        a = cos(theta)^2/2/sigx^2+sin(theta)^2/2/sigy^2;
        b = -sin(2*theta)/4/sigx^2+sin(2*theta)/4/sigy^2;
        c = sin(theta)^2/2/sigx^2+cos(theta)^2/2/sigy^2;
        
        [xdata,ydata] = meshgrid(linspace(-1,1,size(xydata,2)).*size(xydata,2)/2,linspace(-1,1,size(xydata,1)).*size(xydata,1)/2);
        
        fit = A.*exp(-(a.*(xdata-x0).^2+2*b.*(xdata-x0).*(ydata-y0)+c.*(ydata-y0).^2))+B;
% %         
%         subplot(1,2,1)
%         imagesc(xydata);
%         subplot(1,2,2)
%         imagesc(fit);
%         drawnow;
% %         
        sse = sum(sum((fit-xydata).^2));
    end


end
