function [pix,pix_ML] = angle2pix(display,ang)
%pix = angle2pix(display,ang)
%
%converts visual angles in degrees to pixels.
%
%Inputs:
%display.dist (distance from screen (cm))
%display.width (width of screen (cm))
%display.resolution (number of pixels of display in horizontal direction)
%
%ang (visual angle)
%
%Warning: assumes isotropic (square) pixels

%Written 11/1/07 gmb zre

diagonal=display.diag;
resolution=display.res;
diagonal_pix=sqrt(sum(resolution.^2));
width=sqrt(diagonal^2/(1+(resolution(2)/resolution(1))^2));

%Calculate pixel size
pixSize = width/resolution(1);   %cm/pix
sz = 2*display.dist*tan(pi*ang/(2*180));  %cm
pix = sz/pixSize;   %pix 

%Calculate pixel size according to MonkeyLogic
screenangle=(atand((display.diag/2)/display.dist)*2);
ppd_ML = diagonal_pix/screenangle; %pixels per degree according to MonkeyLogi
pix_ML = ppd_ML.*ang;

return