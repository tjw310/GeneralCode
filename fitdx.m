function [est,flag,sse] = fitdx(dx,obj,x,y)

dy = obj.e(2);
t = (360-obj.theta)/180*pi;
dxfit = zeros(size(x,1),1);
op = obj.opticCentre;
model = @fun;
st = [0,0];
options = optimset('MaxfunEvals',2000*length(st),'maxiter',2000*length(st));

for j=1:size(dx,2)
[estimates,~,flag(j)] = fminsearch(model,st,options);
    if flag(j)==1
        sse(j) = model(estimates);
        est(j,:) = estimates;
    else
        sse(j) = NaN;
        est(j,:) = NaN;
    end
end

      function sse = fun(params)
            xa = params(1);
            ya = params(2);

            dxfit = (dy.*(x*cos(t(j))-y*sin(t(j))-op(1)).*(ya*cos(t(j))+xa*sin(t(j))-xa)-xa*cos(t(j))+ya*sin(t(j))+xa)./...
                (1-dy*((x-xa)*sin(t(j))+(y-ya)*cos(t(j))+ya));
            sse = nansum((dxfit-dx(:,j)).^2);
      end

end
        
        