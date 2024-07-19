% add columns to imagematrix
% column luminances are the same as the top pixel of the most lateral column 
% on the side where the new ones are added

function imagematrix=addcols(imagematrix,ntoadd,where)

if strcmpi(where,'left')
pixofinterest=imagematrix(1,1,:);
elseif strcmpi(where,'right')
 pixofinterest=imagematrix(1,end,:);   
else
    error('unclear where to add')
end

newcols=ones(size(imagematrix,1),ntoadd,size(imagematrix,3)).*pixofinterest;

if strcmpi(where,'left')
   newimage=[newcols imagematrix];
elseif strcmpi(where,'right')
   newimage=[imagematrix newcols]; 
else
    error('unclear where to add')
end

imagematrix=newimage;

end