function [theta,phi] = angletestv2(u,v,w)

r = u^2+v^2+w^2

r1 = sqrt(u^2+w^2)

r2 = sqrt(v^2+w^2)

phi = acos(w/r2)/pi*180

phi = atan2(v,w)/pi*180

theta = acos(w/r1)/pi*180

theta = atan2(u,w)/pi*180

end