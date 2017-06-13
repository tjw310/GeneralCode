function [ proj ] = filteringv4(proj,param )
%FILTERING with x,y,z centred on object with u offset on optic axis, v=z. 
[uu,vv] = meshgrid(param.u,param.v);

w = (1+uu.*sqrt(param.R^2-param.D^2)./(param.D^2+vv.^2))./sqrt(param.D^2+uu.^2 + vv.^2);

for i=1:param.nProj
    proj(:,:,i) = proj(:,:,i).*w';
end

if param.parker == 1
	proj = ParkerWeight(proj,param);
end 

filt_len = max(64,2^nextpow2(2*param.nPx));
[ramp_kernel] = ramp_flat(filt_len);

d = 1; % cut off (0~1)
[filt] = Filter(param.filter, ramp_kernel, filt_len, d);
filt = repmat(filt',[1 param.nPx]);

for i=1:param.nProj
    
    fproj = (zeros(filt_len,param.nPx,'single'));
    
    %fproj(filt_len/2-param.nPx/2+1-param.opticCentreP(1):filt_len/2+param.nPx/2+param.opticCentreP(1),:) = proj(:,:,i);
    
    fproj(filt_len/2-param.nPx/2+1:filt_len/2+param.nPx/2,:) = proj(:,:,i);
    
%     subplot(1,2,1)
%     imagesc(fproj);
    
    fproj = fft(fproj);   
    
    fproj = fproj.*filt;
    
    fproj = (real(ifft(fproj)));
    
%     subplot(1,2,2)
%     imagesc(fproj);
%     drawnow;

    
    if param.parker == 1
        proj(:,:,i) = fproj(end/2-param.nPx/2+1+param.opticCentreP(1):end/2+param.nPx/2+param.opticCentreP(1),:)/2/param.du*(2*pi/(180/param.dang))/2;
    else
        %proj(:,:,i) = fproj(end/2-param.nPx/2+1+param.opticCentreP(1):end/2+param.nPx/2+param.opticCentreP(1),:);
        proj(:,:,i) = fproj(end/2-param.nPx/2+1:end/2+param.nPx/2,:);
    end
end


end

function [h, nn] = ramp_flat(n)
nn = [-(n/2):(n/2-1)]';
h = zeros(size(nn),'single');
h(n/2+1) = 1 / 4;
odd = mod(nn,2) == 1;
h(odd) = -1 ./ (pi * nn(odd)).^2;
end


function [filt] = Filter(filter, kernel, order, d)

f_kernel = abs(fft(kernel))*2;
filt = f_kernel(1:order/2+1)';
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

switch lower(filter)
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'  
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        filter
        error('Invalid filter selected.');
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt , filt(end-1:-1:2)];    % Symmetry of the filter
return

end