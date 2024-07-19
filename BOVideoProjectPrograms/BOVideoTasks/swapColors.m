function imagematrix=swapColors(imagematrix)

for chi=1:3
   nonzerolums=unique(imagematrix(:,:,chi)); 
   if numel(nonzerolums)>1
   if numel(nonzerolums)>2, 
       if numel(nonzerolums)==3
       warning('found 3 luminances. Assuming middle one is background grey...'); 
       nonzerolums=[min(nonzerolums) max(nonzerolums)];
       else
          error('can not swap'); 
       end
   end
   
   imtemp=imagematrix(:,:,chi);
   imtemp(imagematrix(:,:,chi)==nonzerolums(1))=nonzerolums(2);
   imtemp(imagematrix(:,:,chi)==nonzerolums(2))=nonzerolums(1);
   imagematrix(:,:,chi)=imtemp;
   end
end

end