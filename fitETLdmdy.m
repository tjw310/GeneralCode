function estimates = fitETLdmdy(images)

x = images.peaksXsub;
z = images.peaksZsub;
y = circshift(x,[0,100])*images.pxSz/images.mtz0;

op1 = images.opticCentre(1);
op2 = images.opticCentre(2);

start_guesses = [images.dmdy];

figure;

model = @fun;
estimates = fminsearch(model,start_guesses);



function sse = fun(params)
    dmdy = params(1);
    if isempty(y)
        y = yGuess;
    end
    
    xTrue = (x-op1).*(1-dmdy.*y)+op1;
    zTrue = (z-op2).*(1-dmdy.*y)+op2;
    %xTrue = xTrue-images.AoRmotion-images.AoRcentreX;
    
    scatter(x,z,'r');hold on; scatter(xTrue,zTrue,[],y); drawnow; hold off; pause(1);
    sse = nanstd(zTrue(:));
    
end


end
    
    
    
    
    