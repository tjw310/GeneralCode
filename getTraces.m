    function [x,z,bool] = getTraces(obj)
        if isempty(obj.peaksX)
            [x,z,~,~] = obj.findBeadsAllSlices;
        else
        x = obj.peaksX; z = obj.peaksZ;
        end

        bool(1:size(x,1),1) = ones(1,size(x,1));
        for i=1:size(x,2)-1
            [xpair,zpair] = findClosestPairs(x(:,i),z(:,i),x(:,i+1),z(:,i+1));

            x(:,i+1) = xpair; z(:,i+1) = zpair;
            
            bool(:,i+1) = ~isnan(x(:,i+1));            
            x(isnan(x(:,i+1)),i+1) = x(isnan(x(:,i+1)),i);
            z(isnan(z(:,i+1)),i+1) = z(isnan(z(:,i+1)),i);
        end
        
        x(bool==0) = NaN;
        z(bool==0) = NaN;
    end