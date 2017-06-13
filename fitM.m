function fitAugust24(xOut,zOut,progX)  


st = [-.05,-900];
model = @fitM;
estimates = fminsearch(model,st);



function fitM(params)
    A = params(1);
    B = params(2);

    fit = A.*sqrt(xOut.^2+zOut.^2)+B;

    scatter(1:length(progX),progX); hold on; plot(fit,'r');hold off; drawnow;

    sse = sum((fit-progX).^2);
end

end