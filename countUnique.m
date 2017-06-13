function [uniq,counts] = countUnique(u)
%u is a 1D array
k=1;
while ~isempty(u)
    if ~isnan(u(1))
    df = abs(u(2:end)-u(1));
    [mn,lc] = min(df);
    if mn==0
        uniq(k) = u(1);
        counts(k) = length(lc)+1;
        idx = 2:length(u);
        for i=1:length(lc)
            idx = idx(idx~=(lc(i)+1));
        end
        u = u(idx);
    else
        uniq(k) = u(1);
        counts(k) = 1;
        u = u(2:end);
    end
    k=k+1;
    else
        u = u(2:end);
    end
end
uniq = uniq.';
counts = counts.';
end
            
    