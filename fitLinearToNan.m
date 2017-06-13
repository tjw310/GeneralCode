function out = fitLinearToNan(array)
%plot(array); drawnow;
% fills in NaN parts with a linear fit to adjacent values
if size(array,1)~=1
    array = array.';
end
%scatter(1:length(array),array); hold on;
out = [];
while ~isempty(array)
    [~,x1] = find(isnan(array),1,'first');
    if isempty(x1)
        out = horzcat(out,array(1,1:end));
        break
    elseif x1==1
        [~,lc] = find(~isnan(array),1,'first');
        array(1,1:lc-1) = array(lc);
    else
    y1 = array(1,x1-1);
    out = horzcat(out,array(1,1:x1-1));    
    array = array(x1:end);
    [~,x2] = find(~isnan(array),1,'first');
    if isempty(x2)
        y2 = out(1,1);
        out = horzcat(out,linspace(out(end),out(end),length(array)));
        array = [];
    else
        y2 = array(1,x2);
        grad = (y2-y1)./(x2+1);
        inter = y2 - grad*x2;
        yfit = grad.*(1:x2-1)+inter;
        out = horzcat(out,yfit);
        array = array(x2:end);
    end
    
    end
end

%plot(out); drawnow;
%hold off;

end