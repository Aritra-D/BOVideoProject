% add rows to imagematrix
% row luminances are the same as the left pixel of the extreme row 
% on the side where the new ones are added

function imagematrix=addrows(imagematrix,ntoadd,where)

if strcmpi(where,'top')
pixofinterest=imagematrix(1,1,:);
elseif strcmpi(where,'bottom')
 pixofinterest=imagematrix(end,1,:);   
else
    error('unclear where to add')
end

newrows=ones(ntoadd,size(imagematrix,2),size(imagematrix,3)).*pixofinterest;

if strcmpi(where,'top')
   newimage=[newrows;imagematrix];
elseif strcmpi(where,'bottom')
   newimage=[imagematrix;newrows]; 
else
    error('unclear where to add')
end

imagematrix=newimage;

end