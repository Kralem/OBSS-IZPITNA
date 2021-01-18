clear all;
%----READ THE INPUT IMAGE----
image = imread ('0005.png');
sigma = min(size(image)) * 0.005;
%----SHOW ORIGINAL IMAGE----
figure, imshow(image);
title('Original image');
%----THRESHOLD AND GAUSS----
T_Low = 0.075;
T_High = 0.175;
gf = imgaussfilt(im2double(image), sigma);
%----SOBEL MASK FOR FILTERING FROM SLIDES----

%Mx = [-1 0 1; -1 0 1; -1 0 1]; %i used the prewitt mask for my first attempt
%My = [-1 -1 -1; 0 0 0; 1 1 1]; %the results weren't good, but they weren't terrible. results were better with sobel 

Mx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
My = [1, 2, 1; 0, 0, 0; -1, -2, -1];
gx = conv2(gf, Mx, 'same');
gy = conv2(gf, My, 'same');
%----ANGLE----
angle = atan2 (gy, gx);
angle = angle*180/pi;
width=size(gf,1);
height=size(gf,2);
%----ADJUST NEGATIVE ANGLES----
for i=1:width
    for j=1:height
        if (angle(i,j)<0) 
            angle(i,j)=360+angle(i,j);
        end
    end
end
angles2=zeros(width, height);
%----ADJUST ANGLES TO NEAREST 0, 45, 90 OR 135 DEGREES----
for i = 1  : width
    for j = 1 : height
        if ((angle(i, j) >= 0 ) && (angle(i, j) < 22.5) || (angle(i, j) >= 157.5) && (angle(i, j) < 202.5) || (angle(i, j) >= 337.5) && (angle(i, j) <= 360))
            angles2(i, j) = 0;
        elseif ((angle(i, j) >= 22.5) && (angle(i, j) < 67.5) || (angle(i, j) >= 202.5) && (angle(i, j) < 247.5))
            angles2(i, j) = 45;
        elseif ((angle(i, j) >= 67.5 && angle(i, j) < 112.5) || (angle(i, j) >= 247.5 && angle(i, j) < 292.5))
            angles2(i, j) = 90;
        elseif ((angle(i, j) >= 112.5 && angle(i, j) < 157.5) || (angle(i, j) >= 292.5 && angle(i, j) < 337.5))
            angles2(i, j) = 135;
        end
    end
end
%----MAGNITUDE----
magnitude = (gx.^2) + (gy.^2);
magnitude = sqrt(magnitude);
thinned = zeros (width, height);
%----EDGE THINNING AKA GRADIENT NON MAXIMUM SUPPRESSION----
for i=2:width-1
    for j=2:height-1
        if (angles2(i,j)==0)
            thinned(i,j) = (magnitude(i,j) == max([magnitude(i,j), magnitude(i,j+1), magnitude(i,j-1)]));
        elseif (angles2(i,j)==45)
            thinned(i,j) = (magnitude(i,j) == max([magnitude(i,j), magnitude(i+1,j-1), magnitude(i-1,j+1)]));
        elseif (angles2(i,j)==90)
            thinned(i,j) = (magnitude(i,j) == max([magnitude(i,j), magnitude(i+1,j), magnitude(i-1,j)]));
        elseif (angles2(i,j)==135)
            thinned(i,j) = (magnitude(i,j) == max([magnitude(i,j), magnitude(i+1,j+1), magnitude(i-1,j-1)]));
        end
    end
end
thinned = thinned.*magnitude;
%----SHOW IMAGE AFTER THINNING----
figure, imshow(thinned);
title('Detected and thinned edges');
%----HYSTERESIS THRESHOLDING----
T_Low = T_Low * max(max(thinned));
T_High = T_High * max(max(thinned));
T_res = zeros (width, height);
for i = 1  : width
    for j = 1 : height
        if (thinned(i, j) < T_Low)
            T_res(i, j) = 0;
        elseif (thinned(i, j) > T_High)
            T_res(i, j) = 1;
        %----8 CONNECTIVITY----
        elseif ( thinned(i+1,j)>T_High || thinned(i-1,j)>T_High || thinned(i,j+1)>T_High || thinned(i,j-1)>T_High || thinned(i-1, j-1)>T_High || thinned(i-1, j+1)>T_High || thinned(i+1, j+1)>T_High || thinned(i+1, j-1)>T_High)
            T_res(i,j) = 1;
        end
    end
end
edge_final = T_res;
%----SHOW FINAL RESULT OF CANNY EDGE----
figure, imshow(edge_final);
title('Final result of canny edge detection');
imwrite(edge_final,'myCanny.png');