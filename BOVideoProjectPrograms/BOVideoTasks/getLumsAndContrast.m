%this function needs at least 2 of min luminance, max luminance and
%lumcontrast, and calculates the remaining parameter. The unknown parameter
%should be NaN.
%Formula used:
%lumcontrast=(maxlum-minlum)/(maxlum+minlum)

function [minlum,maxlum,lumcontrast,bglum]=getLumsAndContrast(minlum,maxlum,lumcontrast,bglum)

if ~isnan(minlum) && ~isnan(maxlum) && ~isnan(lumcontrast)
  return; 
elseif isnan(minlum)+isnan(maxlum)>1 && isnan(bglum)
    error('if minlum and maxlum are unknown, give bglum')
elseif isnan(minlum)+isnan(maxlum)>1
    A=(1-lumcontrast)/(1+lumcontrast);
    minlum=(2*bglum*A)/(1+A);
    maxlum=((1+lumcontrast)*minlum)/(1-lumcontrast);
elseif isnan(minlum)
    minlum=((1-lumcontrast)*maxlum)/(1+lumcontrast);
elseif isnan(maxlum)
    maxlum=((1+lumcontrast)*minlum)/(1-lumcontrast);
else isnan(lumcontrast)
    lumcontrast=(maxlum-minlum)/(maxlum+minlum);
end

end


