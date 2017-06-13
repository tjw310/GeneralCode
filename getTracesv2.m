    function [x,z,bool] = getTracesv2(obj)
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

%         bool(1:size(x,1),1) = ones(1,size(x,1));
%         for i=1:size(x,2)-1
%             [xpair,zpair,exX,exZ] = findClosestPairs(x(:,i),z(:,i),x(:,i+1),z(:,i+1));
%             x(:,i+1) = xpair; z(:,i+1) = zpair;      
%             bool(:,i+1) = ~isnan(x(:,i+1));
%            
%             
%             for l=1:size(x,1)
%                 if isnan(x(l,i+1))
%                     lc = find(~isnan((x(l,2:i)-x(l,1:i-1))),2,'last');
%                     if length(x(l,:))>1 && length(lc)>1    
%                         incX = x(l,lc(2))-x(l,lc(1));
%                         incZ = z(l,lc(2))-z(l,lc(1));
%                         x(l,i+1) = x(l,i)+incX;
%                         z(l,i+1) = z(l,i)+incZ;
%                     else
%                         x(l,i+1) = x(l,i);
%                         z(l,i+1) = z(l,i);
%                     end
%                 end
%                 
%             end
% %             if i==1
% %                 scatter(x(:,i),z(:,i),[],bool(:,i+1)); hold on;
% %             end
% %             scatter(x(:,i+1),z(:,i+1),[],bool(:,i+1)); title(num2str(i)); drawnow;
%           
%             %x(size(x,1)+1:size(x,1)+length(exX),i+1) = exX;
%             %z(size(z,1)+1:size(z,1)+length(exZ),i+1) = exZ;
%         end
%         hold off;
%         
%         x(bool==0) = NaN;
%         z(bool==0) = NaN;
    end