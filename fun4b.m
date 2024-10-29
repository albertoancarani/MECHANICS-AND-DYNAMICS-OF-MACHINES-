function [yp] = fun4b(delta,a,alpha0,alphaprimo)
yp= delta - 2*a*( cos(alpha0)/ cos(alphaprimo) )* (involute(alphaprimo) - involute(alpha0) );
end
