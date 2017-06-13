function [ matrix_out ] = shiftRowsv2( matrix,sortVector )
    for j=1:size(matrix,2)
        if sortVector(j)<=-1
            r = rem(-1*sortVector(j),floor(-1*sortVector(j)));
            matrix_out(:,j) = (1-r).*(circshift(matrix(:,j),[ceil(sortVector(j)),0]))+(r).*(circshift(matrix(:,j),[floor(sortVector(j)),0]));
        elseif sortVector(j)>=1
            r = rem(sortVector(j),floor(sortVector(j)));
            matrix_out(:,j) = (1-r).*(circshift(matrix(:,j),[floor(sortVector(j)),0]))+(r).*(circshift(matrix(:,j),[ceil(sortVector(j)),0]));
        else
            r = abs(sortVector(j));
            if sortVector(j)<0
                matrix_out(:,j) = (1-r).*(circshift(matrix(:,j),[ceil(sortVector(j)),0]))+(r).*(circshift(matrix(:,j),[floor(sortVector(j)),0]));
            else
                matrix_out(:,j) = (1-r).*(circshift(matrix(:,j),[floor(sortVector(j)),0]))+(r).*(circshift(matrix(:,j),[ceil(sortVector(j)),0]));
            end
        end
    end

end

