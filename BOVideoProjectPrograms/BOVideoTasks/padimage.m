% pads image by adding rows and columns evenly on both sides until a total
% of nrows and ncols is attained
function imagematrix=padimage(imagematrix,nrows,ncols)

disp('')

rowstoadd=nrows-size(imagematrix,1);
colstoadd=ncols-size(imagematrix,2);

if colstoadd<0 || rowstoadd<0
   error('image already too large'); 
end
    
if colstoadd>0
tempim=addcols(imagematrix,floor(colstoadd/2),'left');
tempim=addcols(tempim,colstoadd-floor(colstoadd/2),'right');

imagematrix=tempim;
end
if rowstoadd>0
tempim=addrows(tempim,floor(rowstoadd/2),'top');
tempim=addrows(tempim,rowstoadd-floor(rowstoadd/2),'bottom');

imagematrix=tempim;
end


if size(imagematrix,1)~=nrows || size(imagematrix,2)~=ncols
   error('result incorrect'); 
end

end