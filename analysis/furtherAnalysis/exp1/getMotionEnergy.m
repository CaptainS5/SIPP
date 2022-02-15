function netMotionEnergy = getMotionEnergy(imageXYT, tw, screen)
% adapted from the code from http://www.georgemather.com/Model.html;
% the original code was 1D, adapted to 2D but only vertical dimension using
% the spatial and temporal filter from Kiani, Hanks, & Shadlen, 2008.
% Input:
% --imageXYT, each row = one y value, each column = one x value, the third
%       dimension is time in frames; spatial positions are in pixels, each cell
%       is one pixel;
% Output: 
% --net motion energy on the vertical dimension, up is positive.
%
% Explanation from the original code as below:
% Implementation of the Adelson & Bergen (1985) motion energy model.
% Example Matlab code by George Mather, University of Sussex, UK, 2010.
%
% An expanded version of the code is used by George Mather and Kirsten
% Challinor as part of a Wellcome Trust funded research project. Initial
% code was partially based on a tutorial forming part of a short course at
% Cold Spring Harbor Laboratories, USA.
%
% This script is part of an online guide to implementing the Adelson-Bergen
% motion energy model:
%
% http://www.lifesci.sussex.ac.uk/home/George_Mather/EnergyModel.htm
%
% It is free for non-commercial educational use, with appropriate
% acknowledgement.
%
% The script requires a variable called 'stim' to be loaded in
% Step 3b. You can use 'AB15.mat' & 'AB16.mat' or input your own stimulus.
% 
%  
% Adpated by Xiuyun Wu, Jan 27, 2022 (xiuyunwu5@gmail.com)
% units in the image are pixels, calculated from deg

%--------------------------------------------------------------------------
%           STEP 1: Create component spatiotemporal filters
%--------------------------------------------------------------------------

% Define the space axis of the filters
nS=80;              % Number of spatial samples in the filter for both dimensions
[max_x, max_y] = dva2pxl(0.7, 0.7, screen);         % Half-width of filter
% dx = (max_x*2)/nx;  % Spatial sampling interval of filter (deg)

% A row vector holding spatial sampling intervals
x_filt = linspace(-max_x, max_x, nS);
y_filt = fliplr(linspace(-max_y, max_y, nS)); % up is positive
% easier to have the 3D meshfrid before hand with t

% Define the time axis of the filters
max_t = tw/1000; % Duration of impulse response (sec)
nt = round(sec2frm(max_t, screen))+1;         % Number of temporal samples in the filter
% A column vector holding temporal sampling intervals
t_filt = linspace(0, max_t, nt)';

% the 3D meshgrid for the filter
[X, Y, T] =  meshgrid(x_filt, y_filt, t_filt);

% Spatial filter parameters
[sigmaC, ] = dva2pxl(0.35, 0, screen);
[sigmaG, ] = dva2pxl(0.05, 0, screen);
alpha = atan(Y/sigmaC);

% Spatial filter response
gauss = exp(-X.^2/2/sigmaG.^2);          %Gaussian envelope
evenS = cos(alpha).^4 .* cos(4*alpha) .* gauss;   %Even Gabor, corresponds to f1
oddS = cos(alpha).^4 .* sin(4*alpha) .* gauss;    %Odd Gabor, corresponds to f2
% for ii = 1:size(evenS, 3)
%     imshow(evenS(:, :, ii)*255, 'InitialMagnification', 1000)
% end

% Temporal filter parameters
k = 60;    % Scales the response into time units
slow_n = 5;
fast_n = 3;

% Temporal filter response
slow_t = (k*T).^slow_n .* exp(-k*T) .* (1/factorial(slow_n) - (k*T).^2/factorial(slow_n+2));
% corresponds to g2
fast_t = (k*T).^fast_n .* exp(-k*T) .* (1/factorial(fast_n) - (k*T).^2/factorial(fast_n+2));
% corresponds to g1

% Step 1c: combine space and time to create spatiotemporal filters
% element-wise multiplication of X, Y, and T
f1g2 = slow_t .* evenS;    %SE/TS
f1g1 = fast_t .* evenS ;   %SE/TF
f2g2 = slow_t .* oddS ;   %SO/TS
f2g1 = fast_t .* oddS ;   % SO/TF

%--------------------------------------------------------------------------
%         STEP 2: Create spatiotemporally oriented filters
%--------------------------------------------------------------------------
down_1 = f2g1+f1g2;      % L1
down_2 = -f2g2+f1g1;     % L2
up_1 = -f2g1+f1g2;    % R1
up_2 = f2g2+f1g1;     % R2

% % have a look at the filters... they look about right!
% for ii = 1:size(down_1, 3)
%     imshow(up_1(:, :, ii)*255, 'InitialMagnification', 1000)
% end

%--------------------------------------------------------------------------
%         STEP 3: Convolve the filters with a stimulus
%--------------------------------------------------------------------------

% Upward responses
resp_up_1=convn(imageXYT, up_1, 'valid');
resp_up_2=convn(imageXYT, up_2, 'valid');

% Downward responses
resp_down_1=convn(imageXYT, down_1, 'valid');
resp_down_2=convn(imageXYT, down_2, 'valid');

%--------------------------------------------------------------------------
%         STEP 4: Square the filter output
%--------------------------------------------------------------------------
resp_down_1 = resp_down_1.^2;
resp_down_2 = resp_down_2.^2;
resp_up_1 = resp_up_1.^2;
resp_up_2 = resp_up_2.^2;

%--------------------------------------------------------------------------
%         STEP 5: Normalise the filter output
%--------------------------------------------------------------------------
% Calc left and right energy
energy_up= resp_up_1 + resp_up_2;
energy_down= resp_down_1 + resp_down_2;

% Calc total energy
total_energy = sum(sum(energy_up))+sum(sum(energy_down));

% Normalise each directional o/p by total output
RR1 = sum(sum(resp_up_1))./total_energy;
RR2 = sum(sum(resp_up_2))./total_energy;
LR1 = sum(sum(resp_down_1))./total_energy;
LR2 = sum(sum(resp_down_2))./total_energy;
%--------------------------------------------------------------------------
%         STEP 6: Sum the paired filters in each direction
%--------------------------------------------------------------------------
up_Total = RR1+RR2;
down_Total = LR1+LR2;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%         STEP 7: Calculate net energy as the R-L difference
%--------------------------------------------------------------------------
netMotionEnergy = squeeze(up_Total - down_Total);
% reformat into a vector

%--------------------------------------------------------------------------
%         SUPPLEMENTARY CODE: Display summary output and graphics
%--------------------------------------------------------------------------
% % Display motion energy statistic
% fprintf('\n\nNet motion energy = %g\n\n', netMotionEnergy);
% 
% % Plot the stimulus
% figure (1)
% imagesc(stim);
% colormap(gray);
% axis off
% caxis([0 1.0]);
% axis equal
% title('Stimulus');
% 
% % Plot the output:
% %   Generate motion contrast matrix
% energy_opponent = energy_up - energy_down; % L-R difference matrix
% [xv yv] = size(energy_down); % Get the size of the response matrix
% energy_flicker = total_energy/(xv * yv); % A value for average total energy
% 
% % Re-scale (normalize) each pixel in the L-R matrix using average energy.
% motion_contrast = energy_opponent/energy_flicker;
% 
% % Plot, scaling by max L or R value
% mc_max = max(max(motion_contrast));
% mc_min = min(min(motion_contrast));
% if (abs(mc_max) > abs(mc_min))
%     peak = abs(mc_max);
% else
%     peak = abs(mc_min);
% end
% 
% figure (2)
% imagesc(motion_contrast);
% colormap(gray);
% axis off
% caxis([-peak peak]);
% axis equal
% title('Normalised Motion Energy');
%--------------------------------------------------------------------------
end