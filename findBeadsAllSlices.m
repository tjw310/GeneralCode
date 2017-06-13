function [X,Z] = findBeadsAllSlices(obj)
            for i=1:size(obj.projections,3)
                im = obj.projections(:,:,i).';
                [x,z] = FastPeakFind(im,1);
                length(x)
                X(1:length(x),i) = x;
                Z(1:length(z),i) = z;
                disp(i);
                %imagesc(im); hold on; scatter(x,z,'g'); hold off; drawnow;
            end
            X = X-size(obj.projections,1)/2;
            Z = Z-size(obj.projections,2)/2;
end