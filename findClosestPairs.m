function [bxPaired,bzPaired,exX,exZ,mna,mnb,lca,lcb,pairs] = findClosestPairs(ax,az,bx,bz)
%finds closest pairs between points in array a, and array b

for i=1:length(ax)
    ra = sqrt((ax(i)-bx).^2+(az(i)-bz).^2);
    rb = sqrt((ax-bx(i)).^2+(az-bz(i)).^2);
    
    [mna(i),lca(i)] = min(ra);
    [mnb(i),lcb(i)] = min(rb);
end

idx = lcb(lca); %if both have the same minimum pair then would have linear index

bool = idx==(1:length(ax)); %return 1 is have same pair, and 0 if don't

bool2 = mna<25; %2nd test function for min dist<25 pixels

bool = bool.*bool2;

pairs(:,1) = ax(bool==1); bxp = bx(lca); pairs(:,2) = bxp(bool==1);
pairs(:,3) = az(bool==1); bzp = bz(lca); pairs(:,4) = bzp(bool==1);

bxPaired = bxp; exX = bxPaired(bool.'==0); bxPaired(bool.'==0) = NaN;
bzPaired = bzp; exZ = bzPaired(bool.'==0); bzPaired(bool.'==0) = NaN;

end
    

