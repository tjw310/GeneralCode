function [ dmpdy,dmedy ] = getdmedyANDdmpdy(y1,y2,dMp,yp1,yp2,dMe)
%yp1,yp2 are the piezo locations of 180 degree bead images, at fixed ETL
%current and associated change in magnification factor dMp. i.e they are
%the true bead positions.
%ye1, ye2 are the piezo locations for 180 degrees but with different ETL
%currents, etl1, etl2 and associated change in magnification dMe.

%y1,y2,yp1,yp2 are values in air!

n=1.3325;

y1 = n*y1; y2=n*y2; yp1=n*yp1; yp2=n*yp2;

dmpdy = (dMp-1)./(dMp.*y2-y1)

dmedy = (dMe*(1-dmpdy*yp2)-1+dmpdy*yp1)/...
    (dMe*(y2-yp2)*(1-dmpdy*yp2)-(y1-yp1)*(1-dmpdy*yp1))

deltaM = (1-dmpdy*yp1)*(1-dmedy*(y1-yp1))./((1-dmpdy*yp2)*(1-dmedy*(y2-yp2)))

% if rem(dMe,deltaM)~=0
%     error('formula wrong');
% end

deltaY = (yp2-yp1)/2;

dmedy2 = (dMe*(1-dmpdy*deltaY)-1-dmpdy*deltaY)/...
    (dMe*(y2-yp2)*(1-dmpdy*deltaY)-(y1-yp1)*(1+dmpdy*deltaY))

deltaM2 = (1-dmpdy*yp1)*(1-dmedy2*(y1-yp1))./((1-dmpdy*yp2)*(1-dmedy2*(y2-yp2)))

end

