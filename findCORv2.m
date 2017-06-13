function [ optimalShift,model ] = findCORv2( sinogram )
% Written by T.Watson 16/10/2015. Adapted from algorithm presented in
% Reliable Method for calculating centre of rotation in parallel beam
% tomography 2014

shiftguess = 0;

model = @fun;
%options = optimset('MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-20,'TolX',1e-20);
options = optimset('MaxFunEvals',2000,'MaxIter',2000);
[optimalShift] = fminsearch(model,shiftguess,options);


    function [sse] = fun(s)
        %copy sinogram. sinogram is from 0-180 degrees. Flip and shift and add to
        %other half to form full sinogram shifted about horizontal centre.

        shiftedSino = shiftRowsv2(sinogram,repmat(s,[1,size(sinogram,2)]));

        copy = flipud(shiftedSino);

        full = [shiftedSino,copy];

        ft = fftshift(fft2(full));

        dx = 1;

        %max radial extent is half the sinogram x dimension
        r = size(full,1)/2;

        dtheta = 2*pi/size(full,2);

        u = linspace(-1/(2*dx),1/(2*dx),size(full,1));

        n = linspace(-1/(2*dtheta),1/(2*dtheta),size(full,2));


        [n2,u2] = meshgrid(n,u);

        mask = zeros(size(full,1),size(full,2));

        mask(or(n2>abs(r.*u2),n2<-abs(r.*u2)))=1;
        ft = abs(ft).*mask;

        sse = sum(ft(:))/sum(mask(:));
    end

end

