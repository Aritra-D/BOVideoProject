% This function generates a colored square in a background scene with coordinates wrt scene and rotation angle specified
% [x,y] = [0,0 upper left corner of the scene]
function scene = addOrientedSquare(scene,squarewidth,squareheight,x_sqr,y_sqr,rotangle,color,bgval)
close all;

if ~exist('scene','var'); scene = ones(1080,1920).*128; end

% Define circle parameters
scene_size = size(scene); % Size of the output image

if ~exist('squarewidth','var'); squarewidth = 118; end
if ~exist('squareheight','var'); squareheight = 118; end
if ~exist('x_sqr','var'); x_sqr = round(scene_size(2)/2); end % initialize x-coordinate (abscissa) as center of the screen dimension
if ~exist('y_sqr','var'); y_sqr = round(scene_size(1)/2); end % initialize y-coordinate (ordinate) as center of the screen dimension
if ~exist('rotangle','var'); rotangle = 45; end
if ~exist('color','var'); color = 180; end
if ~exist('bgval','var'); bgval = 128; end

cx = round(scene_size(2)/2); % x-coordinate of the circle center
cy = round(scene_size(1)/2); % y-coordinate of the circle center
r = 315;  % Radius of the circle for a 1920 x 1080 display Ref Circle Radius can be set at 315 px.

% % Initialize a blank image (black)
% scene = ones(image_size)*bgval;

% % Starting points for the Midpoint Circle Algorithm
% x = r;
% y = 0;
% decision = 1 - r; % Initial decision parameter
% 
% % Draw the circle outline using the Midpoint Circle Algorithm
% while x >= y
%     % Plot points in all eight octants
%     if cx + x <= scene_size(2) && cy + y <= scene_size(1)
%         scene(cy + y, cx + x) = 0;
%     end
%     if cx - x > 0 && cy + y <= scene_size(1)
%         scene(cy + y, cx - x) = 0;
%     end
%     if cx + x <= scene_size(2) && cy - y > 0
%         scene(cy - y, cx + x) = 0;
%     end
%     if cx - x > 0 && cy - y > 0
%         scene(cy - y, cx - x) = 0;
%     end
%     if cx + y <= scene_size(2) && cy + x <= scene_size(1)
%         scene(cy + x, cx + y) = 0;
%     end
%     if cx - y > 0 && cy + x <= scene_size(1)
%         scene(cy + x, cx - y) = 0;
%     end
%     if cx + y <= scene_size(2) && cy - x > 0
%         scene(cy - x, cx + y) = 0;
%     end
%     if cx - y > 0 && cy - x > 0
%         scene(cy - x, cx - y) = 0;
%     end
% 
%     % Move to the next pixel
%     y = y + 1;
% 
%     % Update the decision parameter to choose between East and South-East pixels
%     if decision <= 0
%         decision = decision + 2 * y + 1;
%     else
%         x = x - 1;
%         decision = decision + 2 * (y - x) + 1;
%     end
% end
% 
% % Draw octant section lines
% for angle = 0:45:315
%     % Convert angle to radians
%     theta = deg2rad(angle);
% 
%     % Calculate the end point of the line based on the circle's radius
%     x_end = round(cx + r * cos(theta));
%     y_end = round(cy - r * sin(theta));
% 
%     % Use Bresenham's line algorithm to draw the line from the center to the edge
%     % Initialize line drawing variables
%     dx = abs(x_end - cx);
%     dy = abs(y_end - cy);
%     sx = sign(x_end - cx);
%     sy = sign(y_end - cy);
%     err = dx - dy;
% 
%     % Bresenham's line drawing loop
%     x = cx;
%     y = cy;
%     while true
%         % Plot the current pixel
%         if x > 0 && x <= scene_size(2) && y > 0 && y <= scene_size(1)
%             scene(y, x) = 0;
%         end
% 
%         % Check if the end of the line is reached
%         if x == x_end && y == y_end
%             break;
%         end
% 
%         % Update the error term and the coordinates
%         e2 = 2 * err;
%         if e2 > -dy
%             err = err - dy;
%             x = x + sx;
%         end
%         if e2 < dx
%             err = err + dx;
%             y = y + sy;
%         end
%     end
%     % figure()
%     % imshow(scene./255);
% end



% Display the result
% imshow(scene./255);


trackRotationPointVal = 254;
scene2 = ones(squarewidth,squareheight)*color;
rotation_point = [1 ceil(squarewidth/2)]; % mid-point of the top side of original square (0 degree) 
scene2(rotation_point(1),rotation_point(2)) = trackRotationPointVal; % make the totation point coordinate (pixel) white to track
rotatedImage = imrotate(scene2,rotangle);
[row_coordinate_trackedRotationPoint,col_coordinate_trackedRotationPoint] = find(rotatedImage==trackRotationPointVal);
rotatedImage(row_coordinate_trackedRotationPoint,col_coordinate_trackedRotationPoint)=color;
% figure(); imshow(scene2./255);
% figure()
% imshow(rotatedImage./255)
% rotatedImage(row_cordinate_trackedRotationPoint,col_cordinate_trackedRotationPoint) = color;
% figure()
% imshow(rotatedImage./255)


% rotAngleList = 0:45:315;
% 
% for i=1:length(rotAngleList)
%     rotatedImage{i} = imrotate(scene2,rotAngleList(i));
%     figure(i)
%     imshow(rotatedImage{i}./255)
% end



[rows,columns] = find(rotatedImage==color);
% rows = sort(rows);
% columns = sort(columns);

% if ~isequal(size(image),size(rotatedImage))
%     rotatedImage(rotatedImage==0) = bgval;
%     [rows,columns] = find(rotatedImage~=bgval);
% else
%     [rows,columns] = find(rotatedImage==color);
% end

x_sqr = x_sqr + round(cx + r * cos(deg2rad(rotangle)));
y_sqr = y_sqr + round(cy - r * sin(deg2rad(rotangle)));

size_image = size(rotatedImage);
image_colcoordinates = round(x_sqr-size(rotatedImage,2)/2):round(x_sqr+size(rotatedImage,2)/2);  %round(x_sqr-size(rotatedImage,2)/2):round(x_sqr+size(rotatedImage,2)/2);
image_rowcoordinates =  y_sqr:y_sqr+size(rotatedImage,1); %round(y_sqr-size(rotatedImage,1)/2):round(y_sqr+size(rotatedImage,1)/2); 

% trim extra rows and columns if necessary to make the image
% dimesnisons fit the grid space

if numel(image_colcoordinates)>size_image(2)
    image_colcoordinates = image_colcoordinates(1:size_image(2));
end

if numel(image_rowcoordinates)>size_image(1)
    image_rowcoordinates = image_rowcoordinates(1:size_image(1));
end

rotationPoint_rowCoordinate = image_rowcoordinates(row_coordinate_trackedRotationPoint);
rotationPoint_colCoordinate = image_colcoordinates(col_coordinate_trackedRotationPoint);

% Translate square so the rotation point is at the origin
translated_image_colcoordinates = x_sqr -rotationPoint_colCoordinate; %repmat(col_cordinate_trackedRotationPoint,[1 size(image_colcoordinates,2)]);
translated_image_rowcoordinates = y_sqr -rotationPoint_rowCoordinate; %repmat(row_cordinate_trackedRotationPoint,[1 size(image_rowcoordinates,2)]);

% Translate the square back to the original position
finalrow_coordinates = image_rowcoordinates +translated_image_rowcoordinates;
finalcolumn_coordinates = image_colcoordinates +translated_image_colcoordinates;



finalrow_coordinates = finalrow_coordinates(finalrow_coordinates>=1 & finalrow_coordinates<=1080);
finalcolumn_coordinates=finalcolumn_coordinates(finalcolumn_coordinates>=1 & finalcolumn_coordinates<=1920);


shape_rowcoordinate_in_image = finalrow_coordinates(rows(rows<=numel(finalrow_coordinates)));
shape_colcoordinate_in_image = finalcolumn_coordinates(columns(columns<=numel(finalcolumn_coordinates)));

% shape_rowcoordinate_in_image = shape_rowcoordinate_in_image(shape_rowcoordinate_in_image>=1);
% shape_colcoordinate_in_image = shape_colcoordinate_in_image(shape_colcoordinate_in_image>=1);


% if isequal(size(shape_rowcoordinate_in_image),size(shape_colcoordinate_in_image))

% scene(shape_rowcoordinate_in_image,shape_colcoordinate_in_image)=color;
% if isequal(size(shape_rowcoordinate_in_image),size(shape_colcoordinate_in_image))
    for iCoordinate=1:size(shape_colcoordinate_in_image,2)
        scene(shape_rowcoordinate_in_image(iCoordinate),shape_colcoordinate_in_image(iCoordinate))=color;
        % scene(rotationPoint_rowCoordinate,rotationPoint_colCoordinate)=trackRotationPointVal;
    end
    scene = scene(1:scene_size(1),1:scene_size(2));
% end

% figure;
% % % % scene(y_sqr,x_sqr) = 0;
% imshow(scene./255)

% end
end
