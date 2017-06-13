function xOut = magscaleX(x,ye,yp,dmpdy,dmedy)

xOut = x./((1-dmpdy.*yp).*(1-dmedy.*(ye-yp)));

end