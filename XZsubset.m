function [Xout,Zout] = XZsubset(x,z,xmn,xmx,zmn,zmx)

for i=1:400
    xa = x(:,i); za = z(:,i);
    xout(i,1:length(xa(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn))))) = xa(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn)));
    zout(i,1:length(za(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn))))) = za(and(and(za<zmx,za>zmn),and(xa<xmx,xa>xmn)));
end

xout(xout==0)=NaN;
zout(zout==0)=NaN;

Xout = nanmean(xout,2).'; Zout = nanmean(zout,2).';

% X2 = Xconc; Z2 = Zconc;
% X2(and(and(Z2<zmx,Z2>zmn),and(X2<xmx,X2>xmn))) = NaN;
% Xout = Xconc(isnan(X2));
% X2 = Xconc;
% Z2(and(and(Z2<zmx,Z2>zmn),and(X2<xmx,X2>xmn))) = NaN;
% Zout = Zconc(isnan(Z2));

figure; scatter(Xout,Zout,[],1:400,'filled'); 
disp(length(Xout));

end
