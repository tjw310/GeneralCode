function s = rotation(mip,shift,angle)

if angle~=0
mip = imrotate(mip,angle,'crop');
end

crSz = floor(size(mip,1)/(2*sqrt(2)));

if shift~=0
    mip = padarray(mip,[500,500],'replicate');
    mip = shiftRowsv2(mip,linspace(shift,shift,size(mip,2)));
    mip = mip(500+1:end-500,500+1:end-500);
end
mip = mip(crSz+1:end-crSz,crSz+1:end-crSz);
%figure; imagesc(mip); drawnow;

sz = 2^8;

fmin = -5; fmax = 5;

g1 = czt_TW2(mip,1,1,sz,fmin/size(mip,1),fmax/size(mip,1));

g2 = czt_TW2(flipud(mip),1,1,sz,fmin/size(mip,1),fmax/size(mip,1));

out = real(g1.*conj(g2));

%out = out./max(out(:));

%figure; imagesc(out); axis square; colorbar(); title(num2str(angle)); drawnow;

q(:,:,1) = out(1:sz/2,1:sz/2);
q(:,:,2) = out(1:sz/2,sz/2+1:end);
q(:,:,3) = out(sz/2+1:end,1:sz/2);
q(:,:,4) = out(sz/2+1:end,sz/2+1:end);

 
% for i=1:4
%     subplot(2,2,i); imagesc(q(:,:,i)); axis square; title(num2str(sum(sum(q(:,:,i)))));
% end
% drawnow;

s(1:4) = sum(sum(q,1),2)./sum(out(:));

end