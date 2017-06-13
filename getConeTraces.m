function beads = getConeTraces( projections )
beads = [];
for i=1:size(projections,3)
    theta = i/size(projections,3)*360;
     nBeads = getBeads(projections(:,:,i).',theta); %new beads in projection image
    for j=1:size(nBeads,2)
        if isempty(beads)
            beads = nBeads(j);
        else
            for k=1:size(beads,2)
                dif(k) = sqrt((beads(k).centreX(end)-nBeads(j).centreX).^2+(beads(k).centreY(end)-nBeads(j).centreY).^2);
            end
            dif(dif==0)=NaN;
            [~,c]=find(dif==min(dif(:)));
            if length(c)>1
                c = max(c);
            end
            vl = dif(c);

            if vl<40
                addPoint(beads(c),nBeads(j));
            elseif i==1
                beads = horzcat(beads,nBeads(j));
            %elseif vl<=40 && vl>=20 && dz>0.04
             %   beads = horzcat(beads,nBeads(j));
            end
        end 
    end
    disp(i/size(projections,3)*100);
    scatterAllPoints(beads,projections(:,:,i).'); drawnow;

end

end

function beads = getBeads(image,theta)
    [x,y]=FastPeakFind(image,1);
    beads = [];         
    if ~isempty(x)
        for i=1:length(x)
            beads = horzcat(beads,beadClass(image(y(i),x(i)),x(i),y(i),[],theta));
        end
    end
end

function scatterAllPoints(beads,image)
imagesc(image);
hold on;
for i=1:size(beads,2)
    x = beads(i).centreX;
    y = beads(i).centreY;
    scatter(x,y);
end
hold off;
end

