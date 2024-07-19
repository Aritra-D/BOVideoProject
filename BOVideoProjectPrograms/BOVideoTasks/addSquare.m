function scene=addSquare(scene,squarewidth,squareheight,x,y,color)
try
scenesize=size(scene);

colcoords=x:x+squarewidth;
rowcoords=y:y+squareheight;
goodis=find(colcoords>=1);
colcoords=colcoords(goodis);
goodis=find(rowcoords>=1);
rowcoords=rowcoords(goodis);

scene(rowcoords,colcoords)=color;

scene=scene(1:scenesize(1),1:scenesize(2));
catch ME
disp('error')
end
end