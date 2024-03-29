function [rdkControl seed] = generateTrialRDKpositions(const, screen, control)
%% generate the position of each dot in all frames
% This is the classical RDK with signal dots + noise dots, controlled by
%   coherence; the current version is Brownian motion
% Each cell is one frame, each row is one dot, position x and y in columns
%   save the seeds for later recovery of the stimuli, output every trial
% *also needs functions from MemToolbox: sd2k, vonmisesrnd
% --rdkControl.dotPos is for plotting the rdk in each frame; each cell is one
%   frame, within each frame each row is one dot, and columns are x and y
%   positions
% --rdkControl.dotDir is the moving directions from the current frame to the next
%   frame, mostly used for debugging; each row is one frame, each column is
%   one dot

% set up RDK
coh = control.rdkCoh;
rdkDuration = control.rdkDuration;
rdkInternalSpeed = control.rdkInternalSpeed;
rdkApertureDir = control.rdkApertureDir;
if rdkApertureDir==0 % moving rightward
    rdkApertureAngle = control.rdkApertureAngle;
else % moving leftward
    rdkApertureAngle = rdkApertureDir-control.rdkApertureAngle;
end
rdkInternalDir = -control.rdkInternalDir; % now the internal direction is fixed within the RDK, irrelative to the aperture direction
% flip for PTB...

% when internal direction is relative to the aperture direction... not in
% use
% % since in PTB it is up-negative, down-positive, and the polar direction is
% % "circular", need to be careful about the correspondence--eventually, an
% % internal dir of 45 means 45 degs above the moving direction. To match
% % this, we need to sort out the values for stimuli display:
% % transfer the relative internal motion direction to absolute direction for
% % display
% if rdkApertureDir==0 % moving rightward
%     rdkInternalDir = control.rdkApertureDir-control.rdkInternalDir; % flip since for PTB it is up-negative, down-positive
% else
%     rdkInternalDir = control.rdkApertureDir+control.rdkInternalDir;
% end

% if const.startExp==1 || const.startExp==0 % actual experiment, translating aperture
    % initialize horizontal aperture movement per frame
    [moveDistanceAperture, ] = dva2pxl(const.rdk.apertureSpeed, const.rdk.apertureSpeed, screen)*screen.refreshRate; % pixel per frame
    % initialize the starting center position relative to center of screen
    [apertureStartDisX, apertureStartDisY] = dva2pxl(rdkDuration*const.rdk.apertureSpeed/2*cos(rdkApertureAngle/180*pi), ...
        rdkDuration*const.rdk.apertureSpeed/2*sin(rdkApertureAngle/180*pi), screen); % the distance between the starting center point and center of screen
    % randomize the trajectory center to randomize the starting position
    [jitterDisX, jitterDisY] = dva2pxl(const.startingPositionJitter, const.startingPositionJitter, screen);
    rdkControl.randCenterX = round(screen.x_mid+(2*rand-1)*jitterDisX);
    rdkControl.randCenterY = round(screen.y_mid-(2*rand-1)*jitterDisY);
    rdkControl.apertureCenterPos{1} = [rdkControl.randCenterX-apertureStartDisX rdkControl.randCenterY+apertureStartDisY]; % the initial starting position, depends on moving direction and speed
% elseif const.startExp==-1 % baseline, static aperture
%     rdkControl.apertureCenterPos{1} = screen.center;
% end

% Define the texture center and texture window for both the whole aperture texture and RDK dot field 
% They should always overlap, regardless whether the aperture is the same
% size of the whole texture (dots move together with the aperture) or the
% aperture is a smaller "hole" translating across the dot field.
% Dots will be plotted relative to the center position of the whole  aperture texture, 
% which may not be the aperture center if the aperture is translating across the dot field
[dotFieldRadiusX, dotFieldRadiusY] = dva2pxl(const.rdk.dotFieldRadius, const.rdk.dotFieldRadius, screen);
% if const.apertureType==1 % aperture translates across the dot field
%     rdkControl.textureCenterPos = screen.center;
%     rdkControl.textureWindow = [rdkControl.textureCenterPos(1)-dotFieldRadiusX, ...
%         rdkControl.textureCenterPos(2)-dotFieldRadiusY, rdkControl.textureCenterPos(1)+dotFieldRadiusX, ...
%         rdkControl.textureCenterPos(2)+dotFieldRadiusY]; % the window to draw aperture texture in
% else % dots move together with the aperture
    rdkControl.textureCenterPos{1} = rdkControl.apertureCenterPos{1};
    rdkControl.textureWindow{1} = [rdkControl.textureCenterPos{1}(1)-dotFieldRadiusX, ...
        rdkControl.textureCenterPos{1}(2)-dotFieldRadiusY, rdkControl.textureCenterPos{1}(1)+dotFieldRadiusX, ...
        rdkControl.textureCenterPos{1}(2)+dotFieldRadiusY]; % the window to draw aperture texture in
% end

% initialize timing parameters
rdkFrames = ceil(sec2frm(rdkDuration, screen));
rdkLifeTime = round(sec2frm(const.rdk.lifeTime, screen));

% initialize RDK dot parameters
[dots.radiusX, ] = dva2pxl(const.rdk.dotRadius, const.rdk.dotRadius, screen);
dots.diameterX = dots.radiusX*2;
% since in dva2pxl we already equate cm and pixels, later they are just used interchageablly, conceptually mixed...

% should already be shuffled at the begining of each experiment, now record
% the exact seed state for each trial
seed = rng;

% Postion dots in a circular aperture using distanceToCenter and
% positionAxis (calculated from theta)
% 1. The initial position (center of the dots) are determined by distance to center and theta, the angle of the distannce vector
dots.distanceToCenterX{1} = dotFieldRadiusX * sqrt((rand(const.rdk.dotNumber, 1))); % distance of dots from center of the aperture
% just use the aperture radius to make sure dots are within the aperture...
% not sure how to use drawdots with transparent aperture; need to create texture for this
theta = 2 * pi * rand(const.rdk.dotNumber,1); % values between 0 and 2pi (2pi ~ 6.28)
rdkControl.dotPos{1} = [cos(theta) sin(theta)] .* [dots.distanceToCenterX{1} dots.distanceToCenterX{1}*screen.pixelRatioWidthPerHeight];

% initialize move distance based on dot speed
[moveDistanceDot, ] = dva2pxl(rdkInternalSpeed, rdkInternalSpeed, screen)*screen.refreshRate; % pixel per frame
moveDistanceDot = repmat(moveDistanceDot, const.rdk.dotNumber, 1);

% initialize dot life time and label time
dots.showTime{1} = rand(1, const.rdk.dotNumber).*ones(1, const.rdk.dotNumber)*rdkLifeTime; % in frames
dots.labelTime(1) = round(sec2frm(const.rdk.labelUpdateTime, screen)); % in frames

% assign dot labels
targetDotsN = round(coh*const.rdk.dotNumber); % number of dots should be moving coherently
dots.label{1} = [ones(targetDotsN, 1); zeros(const.rdk.dotNumber-targetDotsN, 1)]; % target = 1, noise = 0

% initialize moving directions； 0/2pi equals to the horizontal right, rotates CW
moveTheta = 2 * pi * rand(const.rdk.dotNumber, 1); % all random directions except 0/2pi, or the horizontal right
moveTheta(1:targetDotsN, 1) = 0; % assign the signal dots to be coherently moving rightwards
moveTheta = moveTheta+rdkInternalDir/180*pi; % rotate the signal direction to the defined direction
rdkControl.dotDir(1, :) = -moveTheta; % in radians, for this up is positive, down is negative
% movement of each dot from frame N to frame N+1, now initialize for frame1
dots.movement{1} = [cos(moveTheta) sin(moveTheta)].*[moveDistanceDot moveDistanceDot*screen.pixelRatioWidthPerHeight];

% transparent motion noise, fixed label for target and noise dots; 
%   noise dots moving in a new random direction after reappearance,
%   while target dots always have the same moveTheta
% to use Brownian motion, dots.label is updated later in each frame

% generate the dot matrices of the RDK for the whole trial
for frameN = 1:rdkFrames-1
    % update the center position of the translating aperture
%     if const.startExp==1 || const.startExp==0
        rdkControl.apertureCenterPos{frameN+1} = [rdkControl.apertureCenterPos{frameN}(1)+moveDistanceAperture*cos(rdkApertureAngle/180*pi), ...
            rdkControl.apertureCenterPos{frameN}(2)-moveDistanceAperture*sin(rdkApertureAngle/180*pi)];
        rdkControl.textureCenterPos{frameN+1} = rdkControl.apertureCenterPos{frameN+1};
        rdkControl.textureWindow{frameN+1} = [rdkControl.textureCenterPos{frameN+1}(1)-dotFieldRadiusX, ...
            rdkControl.textureCenterPos{frameN+1}(2)-dotFieldRadiusY, rdkControl.textureCenterPos{frameN+1}(1)+dotFieldRadiusX, ...
            rdkControl.textureCenterPos{frameN+1}(2)+dotFieldRadiusY]; % the window to draw aperture texture in
%     elseif const.startExp==-1 % baseline, static aperture
%         rdkControl.apertureCenterPos{frameN+1} = rdkControl.apertureCenterPos{frameN};
%         rdkControl.textureCenterPos{frameN+1} = rdkControl.apertureCenterPos{frameN+1};
%         rdkControl.textureWindow{frameN+1} = rdkControl.textureWindow{frameN};
%     end
    
    if coh==0 % static pattern; if want to use noise pattern, just comment out this "if"
        % no need to update dot position, just copy...
        rdkControl.dotPos{frameN+1} = rdkControl.dotPos{frameN};
        dots.showTime{frameN+1} = dots.showTime{frameN}-1;
    else
        % update dot position
        rdkControl.dotPos{frameN+1} = rdkControl.dotPos{frameN} + dots.movement{frameN}; % without aperture movement
        % update lifetime and labels
        dots.showTime{frameN+1} = dots.showTime{frameN}-1;
        dots.labelTime(frameN+1) = dots.labelTime(frameN)-1;
        
        % initialize for use in the next loop, from frame N+1 to N+2
        dots.movement{frameN+1} = dots.movement{frameN};
        dots.label{frameN+1} = dots.label{frameN};
        rdkControl.dotDir(frameN+1, :) = rdkControl.dotDir(frameN, :);
        
        % renew labels
        if dots.labelTime(frameN+1)<=0
            dots.labelTime(frameN+1)=round(sec2frm(const.rdk.labelUpdateTime, screen));
            % generate new labels
            labelOrder = randperm(size(dots.label{frameN}, 1));
            dots.label{frameN+1}(:, 1) = dots.label{frameN}(labelOrder, 1); % randomly assign new labels
            % update directions
            moveTheta = 2 * pi * rand(const.rdk.dotNumber, 1);
            moveTheta(dots.label{frameN+1}==1, :) = 0; % signal dots moving horizontally to the right
            moveTheta = moveTheta+rdkInternalDir/180*pi; % rotate the signal direction to the defined direction
            rdkControl.dotDir(frameN+1, :) = -moveTheta; % in radians, for this up is positive, down is negative
            % movement of each dot from frame N+1 to frame N+2
            dots.movement{frameN+1} = [cos(moveTheta) sin(moveTheta)].*[moveDistanceDot moveDistanceDot*screen.pixelRatioWidthPerHeight];
        end
        
        % still needs to replace expired dots and move dots out of the aperture into the aperture again, from the opposite edge
%         % 1. Replace dots with expired lifetime
%         expiredDots = find(dots.showTime{frameN+1}' <= 0);
%         if expiredDots
%             dotsN = length(expiredDots);
%             dis2CenterX = dotFieldRadiusX * sqrt((rand(dotsN,1)));
%             theta = 2 * pi * rand(dotsN,1);
%             % generate new positions and update lifetime
%             rdkControl.dotDir(frameN+1, expiredDots) = -theta; % in radians, for this up is positive and down is negative
%             rdkControl.dotPos{frameN+1}(expiredDots, :) = [cos(theta) sin(theta)] .* [dis2CenterX dis2CenterX*screen.pixelRatioWidthPerHeight];
%             dots.showTime{frameN+1}(expiredDots) = rdkLifeTime;
%         end
        % 2. Relocate dots out of the aperture
        dotDist = rdkControl.dotPos{frameN+1}(:, 1).^2 + ...
            ((rdkControl.dotPos{frameN+1}(:, 2)/screen.pixelRatioWidthPerHeight)).^2;
        outDots = find(dotDist>dotFieldRadiusX^2); % all dots out of the aperture
        %     % simply generate new random positions for these dots--to avoid
        %     % opposite direction at the edge of a small aperture (for example,
        %     % signal direction is up, a dot is at the left edge of a 1deg-RDK, two
        %     % steps it will move out of the aperture, then it could look like moving
        %     % downward constantly...)
        %     if outDots
        %         dotsN = length(outDots);
        %         dis2CenterX = dotFieldRadiusX * sqrt((rand(dotsN,1)));
        %         theta = 2 * pi * rand(dotsN,1);
        %         % generate new positions
        %         rdkControl.dotDir(frameN+1, outDots) = -theta; % in radians, for this up is positive and down is negative
        %         rdkControl.dotPos{frameN+1}(outDots, :) = [cos(theta) sin(theta)] .* [dis2CenterX dis2CenterX*screen.pixelRatioWidthPerHeight];
        %         dots.showTime{frameN+1}(outDots) = rdkLifeTime;
        %     end
        %
        % move dots in the aperture from the opposite edge, continue the assigned motion
        rdkControl.dotPos{frameN+1}(outDots, :) = -rdkControl.dotPos{frameN+1}(outDots, :)+dots.movement{frameN}(outDots, :);
    end
    % Replace dots with expired lifetime
    expiredDots = find(dots.showTime{frameN+1}' <= 0);
    if expiredDots
        dotsN = length(expiredDots);
        dis2CenterX = dotFieldRadiusX * sqrt((rand(dotsN,1)));
        theta = 2 * pi * rand(dotsN,1);
        % generate new positions and update lifetime
        rdkControl.dotDir(frameN+1, expiredDots) = -theta; % in radians, for this up is positive and down is negative
        rdkControl.dotPos{frameN+1}(expiredDots, :) = [cos(theta) sin(theta)] .* [dis2CenterX dis2CenterX*screen.pixelRatioWidthPerHeight];
        dots.showTime{frameN+1}(expiredDots) = rdkLifeTime;
    end
end