function [X,Z,Xconc,Zconc] = findBeadsAllSlicesv2(projections)
        Xconc = []; Zconc = [];
        %thres = nanmean(obj.projections(:))-sqrt(nanstd(obj.projections(:)));
        %thres = nanmean(projections(:))-nanstd(projections(:))/5;
        for i=1:size(projections,3)
            im = projections(:,:,i);
            [z,x] = FastPeakFind(im,1);
            X(1:length(x),i) = x;
            Z(1:length(z),i) = z;
            Xconc = horzcat(Xconc,x);
            Zconc = horzcat(Zconc,z);
            disp(i);
            imagesc(im); hold on; scatter(z,x,'g'); hold off; drawnow;
        end
        X = X-size(projections,2)/2;
        Z = Z-size(projections,1)/2;
        Xconc = Xconc -size(projections,2)/2;
        Zconc = Zconc -size(projections,1)/2;
    end  