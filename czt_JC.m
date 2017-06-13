function g = czt_JC(x, k, w, a)
%   James Clegg Edit for phase tilt correction. (Accomplishes same effect
%   as fftshift in ordinary MATLAB fourier transforms). This addtional
%   phase shift ensures that the real and imaginary components of the field
%   overlap correctly.

g = czt(x, k, w, a);
%% From original CZT.m file
[m, n] = size(x); oldm = m;
if m == 1, x = x(:); [m, n] = size(x); end 

%% MN and JHC correction
pcorr = w^((-m+1)/2).^((1:k)')*(a^(k/4));
% disp(size(repmat(pcorr,1,n)));
% disp(size(g))
g = g.*repmat(pcorr,1,n);

%% From original CZT.m file 
if oldm == 1, g = g.'; end

