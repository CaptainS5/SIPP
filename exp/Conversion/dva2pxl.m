function [lengthPixelX lengthPixelY] = dva2pxl(lengthX, lengthY, screen)
% translate length in dva to in pixels
% the length is from fixation to one end, not centered at the fixation
% if want to calculate a length centered at fixation, use length/2 for
% pixels, then multiply 2 to the pixel number calculated from length/2

% but actually it doesn't matter that much for our set-up...
% % first calculate how much cm it is for the visual degree
% lengthXcm = tan(lengthX/180*pi)*screen.dist;
% lengthYcm = tan(lengthY/180*pi)*screen.dist;
% 
% % then calculate how many pixels
% lengthPixelX = round(lengthXcm*screen.ppcX);
% lengthPixelY = round(lengthYcm*screen.ppcY);

% so could just use the ppd
lengthPixelX = round(lengthX*screen.ppd);
lengthPixelY = round(lengthY*screen.ppd);
end