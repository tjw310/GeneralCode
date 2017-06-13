function [u,v,w] = getAxDrVc(alpha,beta)
%gets axes irection vector based on alpha (rotation about y) and beta (tilt
%along z)

alpha = alpha/180*pi;
beta = beta/180*pi;

r1 = sqrt((1-sin(alpha)^2)/(1-sin(alpha)^2*sin(beta)^2));

r2 = sqrt((1-sin(beta)^2)/(1-sin(alpha)^2*sin(beta)^2));

u = r2*sin(alpha);

v = r1*sin(beta);

w = r1*cos(beta);

%scatter3(u,v,w); drawnow;

r = r1^2+r2^2-w^2; %this is a check to get a unit vector

end