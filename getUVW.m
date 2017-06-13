function [u,v,w] = getUVW(theta,phi)

theta = theta/180*pi;
phi = phi/180*pi;

w = (1+(tan(theta)^2)+(tan(phi)^2))^-0.5;

u = w*tan(theta);

v = w*tan(phi);

if or(and(phi<pi/4,phi>0),phi<-pi/4)
        w = -w;
end

r = u^2+v^2+w^2;

end