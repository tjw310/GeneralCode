function fitDiscwithCont(contobj,discobj)

start_guesses = contobj.AoRangle;
model = @fun;
options = optimset('MaxFunEvals',5000*length(start_guesses),'MaxIter',5000*length(start_guesses));
[estimates,~,flag] = fminsearch(model,start_guesses,options);


function sse = fun(params)
        theta = params(1)
        
        dx = discobj.stmotX/contobj.pxSz.*contobj.mtz0*cos(theta/180*pi)-discobj.stmotZ/contobj.pxSz.*contobj.mtz0*sin(theta/180*pi);
        dz = discobj.stmotZ/contobj.pxSz.*contobj.mtz0*cos(theta/180*pi)+discobj.stmotX/contobj.pxSz.*contobj.mtz0*sin(theta/180*pi);
        
        x2 = discobj.peaksXsub+dx;
        z2 = discobj.peaksZsub+dz;
        
        scatter(contobj.peaksXsub,contobj.peaksZsub); hold on;
        scatter(x2,z2); hold off; drawnow;
        
        sse = nansum(sqrt((x2-contobj.peaksXsub).^2+(z2-contobj.peaksZsub).^2));
        
       % sse = nansum((x2-contobj.peaksXsub).^2);
        
       % sse = nansum((z2-contobj.peaksZsub).^2);
        
end

end