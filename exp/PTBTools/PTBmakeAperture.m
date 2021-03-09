function aperture = PTBmakeAperture(const, screen)

Screen('BlendFunction', screen.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[apertureRadiusX, apertureRadiusY] = dva2pxl(const.rdk.apertureRadius, const.rdk.apertureRadius, screen);
imgApt(:, :, 1) = ones(2*apertureRadiusY+1, 2*apertureRadiusX+1)*screen.background; % background layer
transApt = zeros(2*apertureRadiusY+1, 2*apertureRadiusX+1); % transparency layer, now all transparent
[x,y]=meshgrid(-apertureRadiusX:apertureRadiusX, -apertureRadiusY:apertureRadiusY);
x = x/apertureRadiusX;
y = y/apertureRadiusY;
radiusA = sqrt(x.^2+y.^2);
transApt(radiusA>1) = 255; % opaque out of the circle
imgApt(:, :, 2) = transApt;
aperture = Screen('MakeTexture', screen.window, imgApt);

end