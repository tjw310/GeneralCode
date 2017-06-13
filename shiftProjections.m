function projections = shiftProjections(projections,param)
%shift projections relative to optic centre. (so that optic centre is in
%the centre

nr = size(projections,1);
nc = size(projections,2);

c = param.opticCentre(1)./param.pxSz*param.mtz0
r = param.opticCentre(2)./param.pxSz*param.mtz0

if r>0
    r1 = floor(r); r2 = rem(r,1);
else
    r1 = ceil(r); r2 = -1*rem(r,1);
end

if c>0      
    c1 = floor(c); c2 = rem(c,1);
else
    c1 = ceil(c); c2 = -1*rem(c,1);
end
        

projections = padarray(projections,[abs(r1)+1,abs(c1)+1,0],'both');

p1 = circshift(projections,[r1,c1,0]);
p2 = circshift(projections,[r1+1,c1+1,0]);

projections = sqrt(r2^2+c2^2).*p1 + (1-sqrt(r2^2+c2^2)).*p2;

projections = circshift(projections,[-abs(r1)-1,-abs(c1)-1,0]);

projections = projections(1:nr,1:nc,:);


end

