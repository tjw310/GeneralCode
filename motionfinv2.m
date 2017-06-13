function [avDf,df1,df2] = motionfinv2(x1,x2,obj)

ymot = 1.3325*circshift((obj.piezomotion-mean(obj.piezomotion))/obj.pxSz*obj.mtz0,[0,obj.nProj/4]);

model = @fun;
st_g1 = fitSinusoid(x1);
st_g2 = fitSinusoid(x2);

figure;
estimates = fminsearch(model,[st_g1,st_g2])
[~,fitA,fitB,avDf,df1,df2] = model(estimates);

figure; plot(fitA,'b'); hold on; plot(fitB,'g'); plot(ymot,'r'); hold off;
figure; plot(fitA-ymot,'b'); hold on; plot(fitB-ymot,'g'); hold off;

figure; plot(circshift(fitA-ymot,[0,-obj.nProj/4]),'b');
hold on; plot(circshift(fitB-ymot,[0,-obj.nProj/4]),'g'); hold off;

    function [sse,fit1,fit2,avDf,df1,df2] = fun(params)
        a1 = params(1);
        b1 = params(2);
        c1 = params(3);
        a2 = params(4);
        b2 = params(5);
        c2 = params(6);
        
        fit1 = a1*sin(obj.theta/180*pi+b1)+c1;
        fit2 = a2*sin(obj.theta/180*pi+b2)+c2;
        
        df1 = fit1-x1; df2 = fit2-x2;
        
        %df1(abs(df1)>30)=NaN;  df2(abs(df2)>30)=NaN; 
        
        plot(df1); hold on; plot(df2); hold off; drawnow;
        
        avDf = (df1+df2)/2;
        
        sse = nansum(df1.^2)+nansum(df2.^2);
    end
        

end
        