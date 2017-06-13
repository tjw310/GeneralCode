function [xs,zs] = scaleCoords(x,y,z,op,dmdy)
    r =sqrt((x-op(1)).^2+(z-op(2)).^2);
    theta = atan2((z-op(2)),(x-op(1)));
    dm = 1./(1-dmdy.*y);
    dr = r.*(dm-1);
    xs = x+dr.*cos(theta);
    zs = z+dr.*sin(theta);
end