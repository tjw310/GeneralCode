function [rotAngle,shiftAm] = findRotShift(mip,varargin)

if nargin<2
k=1; angles = -5:0.5:5; shifts = -100:25:100;

for i=shifts
    l=1;
    for j=angles
        s(k,l,:) = rotation(mip,i,j);
        l=l+1;
    end
    k=k+1;
    disp(i);
end; 


figure;
for i=1:4
    subplot(2,2,i);
    imagesc(angles,shifts,s(:,:,i));
end

figure; imagesc(angles,shifts,std(s,0,3));

dif = std(s,[],3);
[r,c] = find(dif==min(dif(:)));
rotAngle = angles(c); shiftAm = shifts(r);
sprintf('Initial guess at shift=%f,rot=%f,sse=%f',rotAngle,shiftAm,min(dif(:)))
else
     shiftAm = varargin{1}; rotAngle = varargin{2};
end

start_guesses = [shiftAm,rotAngle];
model = @fun;
est = fminsearch(model,start_guesses);

rotAngle = est(2); shiftAm = -est(1);

    function sse = fun(params)
        A = params(1);
        B = params(2);
        if or(or(A<-100,A>100),or(B<-5,B>5))
            sse=1;
        else
        out = rotation(mip,A,B);
        sse = std(out,0,2);
        sprintf('shift=%f,rot=%f,sse=%f',A,B,sse)
        end
    end
end
