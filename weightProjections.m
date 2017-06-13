function p = weightProjections(obj)

u = obj.x-obj.opticCentre(1)/obj.SubSampFc;

OaDl = obj.nPx/2-abs(obj.AoRcentreX/obj.SubSampFc);
OsOa = (obj.AoRcentreX-obj.opticCentre(1))/obj.SubSampFc;
OsDl = OsOa+OaDl;
D = obj.D;

%R = sqrt(OsOa^2+D^2);
%prefact = R^2./sqrt(R^2+u.^2);
%prefact = (obj.D-u./obj.D.*OsOa)./sqrt(obj.D.^2+u.^2);

proj = squeeze(obj.projections(:,size(obj.projections,1)/2,:));
%pM = proj.*repmat(prefact.',1,size(obj.projections,3));
pM = proj;

gamma = atan2(u,D);
theta = atan2(OsDl,D);
rho = gamma - atan2(OsOa,D);
ep = atan2(OsDl,D)-atan2(OsOa,D);

weight = 0.5*(sin(pi*rho./(-2*ep))+1);

weight(rho<=-ep)=1;

p = pM.*repmat(weight.',1,size(obj.projections,3));

subplot(2,2,1); imagesc(proj);
subplot(2,2,2); imagesc(pM);
subplot(2,2,3); imagesc(p);

end