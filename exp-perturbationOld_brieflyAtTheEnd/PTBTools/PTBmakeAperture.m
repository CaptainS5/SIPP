function aperture = PTBmakeAperture(const, screen, apertureCenterPos)
% apertureCenterPos is in pixel coordinates on the whole screen for apertureType 1 (aperture translates across the dot field), 
% and 0 for apertureType 0 (dots move together with the aperture)

Screen('BlendFunction', screen.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% the whole aperture texture should always cover the dot texture; in normal
% condition they would be the same size, but for a translating aperture on 
% the whole dot field, the aperture will only reveal a small part of the dot field
[apertureRadiusX, apertureRadiusY] = dva2pxl(const.rdk.apertureRadius, const.rdk.apertureRadius, screen);
[dotFieldRadiusX, dotFieldRadiusY] = dva2pxl(const.rdk.dotFieldRadius, const.rdk.dotFieldRadius, screen);

% first, form the background layer of the whole texture, which should be the 
% same size as the dot field
imgApt(:, :, 1) = ones(2*dotFieldRadiusY+1, 2*dotFieldRadiusX+1)*screen.background; % background layer
% +1 because there will be a zero in between...

% initialize transparent layer, now all opaque
transApt = ones(2*dotFieldRadiusY+1, 2*dotFieldRadiusX+1)*255; % transparency layer

% create the meshgrid for the aperture
[x,y]=meshgrid(-apertureRadiusX:apertureRadiusX, -apertureRadiusY:apertureRadiusY);
x = x/apertureRadiusX;
y = y/apertureRadiusY;
radiusA = sqrt(x.^2+y.^2);

% if any(apertureCenterPos~=0) % aperture translates across the dot field
%     % first, create the aperture related values in the transparent layer,
%     % the same size as radiusA, kind of a patch that we will later put into
%     % the transparent layer
%     patchApt = ones(size(radiusA))*255;
%     patchApt(radiusA<=1) = 0; % only within the aperture is transparent
%     
%     patchCenterPos = apertureCenterPos-screen.center+[dotFieldRadiusX, dotFieldRadiusY]+1; % align the position with the current texture center
%     % now the apertureCenterPos is the location in imgApt, row and column idx
%     
%     % now put patchApt to the correct position within the whole texture
%     transApt(patchCenterPos(1)-apertureRadiusY:patchCenterPos(1)+apertureRadiusY, ...
%         patchCenterPos(2)-apertureRadiusX:patchCenterPos(2)+apertureRadiusX) = patchApt;
% else % dots move together with the aperture, the aperture is as large as the dot field
    transApt(radiusA<=1) = 0; % transparent inside the circle
% end
imgApt(:, :, 2) = transApt;
aperture = Screen('MakeTexture', screen.window, imgApt);

end