function motionfin(obj)

x = obj.peaksXsub;
ymot = 1.3325*circshift((obj.piezomotion-mean(obj.piezomotion))/obj.pxSz*obj.mtz0,[0,obj.nProj/4]);

model=@fun; model2=@refine;
st = [0.1,obj.AoRcentreX];
e = fminsearch(model,st); 
e = fminsearch(model2,e)
[~,df,fit] = model2(e);

figure; plot(df);
figure; plot(x); hold on; plot(fit); hold off;

function [sse,df,fit] = fun(p)
        a = p(1);
        b = p(2);
        
         fit = (x-b)*a;
        
        df = fit-ymot;   
        sse = nansum(df.^2);
end

function [sse,df,fit] = refine(p)
        a = p(1);
        b = p(2);
        
        fit = (x-b)*a;
        
        df = fit-ymot;
        df(abs(df)>50)=NaN;
        sse = nansum(df.^2);
end

end
        