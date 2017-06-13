function [ xOut,yOut ] = findClosestLineIntercept( gradients,intercepts )
% searches for point in (x,y) in which it results in the minimum distance
% between a series of crossed lines. i.e if all were perfect it would just
% be at the interception of the lines

x=1:2560;
y = repmat(gradients.',1,length(x)).*repmat(x,length(gradients),1)+repmat(intercepts.',1,length(x));

start_guesses = [x(round(length(x)/2)),mean(y(:))];
model = @fun;
options = optimset('MaxFunEvals',2000*length(start_guesses),'MaxIter',2000*length(start_guesses));
[estimates] = fminsearch(model,start_guesses,options);

xOut = estimates(1);
yOut = estimates(2);

    function [sse] = fun(params)
        xCent = params(1);
        yCent = params(2);
        
        for i=1:size(y,1)
            k = y(i,:);
            r = sqrt((x-xCent).^2+(k-yCent).^2);
            mn(i) = min(r);
        end
        sse = nansum(mn);
    end


end
