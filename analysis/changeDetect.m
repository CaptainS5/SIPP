% FUNCTION to find pursuit onset
% find the change in y(x), by fitting a piecewise linear model
% with one break. x is time (frames), y is the 2D velocity vector
%         / ry      
%   ly\  / 
%      \/
%       (cx,cy)
% It works for the curve with either shape (pointing down or up).
% If you want you can also use this function to find the steady-state 
% phase onset (for that before the break point should only include 
% the open-loop phase, since this only assumes one break point), 
% although a fixed open-loop duration of 140ms works just fine under 
% most conditions. You can also modify the function to fit a piecewise
% linear model with two breaks for finding both onsets.
% History
% this function was written by Dinesh Pai some time ago
% 26-04-2018    XW added some comments
% for questions email xiuyunwu5@gmail.com

function [cx,cy,ly,ry] = changeDetect(x,y)
warning off; 
% initialize parameters, p0, to reasonable values
w = ceil(length(x)/10);% small window
ly = mean(y(1:w)); % y value of the most left point
ry = mean(y(end-w+1:end)); % y value of the most right point
cy = mean(y(w+1:end-w)); % y value of the break point
cx = mean(x(w+1:end-w)); % the time point of the break, which is pursuit onset
p0 = [cx,cy,ly,ry];
% minimize residual
options = optimset('Display', 'off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The original code, which uses the L2 norm, being more sensitive to variations
%   would work well for normal pursuit traces, but maybe not so much for
%   some specific conditions; choose depending on your data
p = lsqnonlin(@myfun,p0,[],[],options); % uncomment this line and comment section below to use
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Instead, using the L1 norm will be less sensitive to changes during pursuit
% % initiation, the smooth L1 norm version defined as below, uncomment to use 
% p = fminunc(@smoothL1, p0, options);
% % but it actually doesn't really matter in my case... so could ignore this
% % part for now; or try further restricting the smooth part... Huber loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cx = p(1); cy = p(2); ly = p(3); ry = p(4);
return;
% inner function for computing residual errors
    function residual = myfun (p)
        cx = p(1); cy = p(2); ly = p(3); ry = p(4);
        residual = zeros(size(x));
        for ii = 1:length(x)
            residual(ii) = y(ii) - evalPWL(x(ii),x(1),ly,cx,cy,x(end),ry);
        end
    end
% smooth L1 norm cost function, not really necessary now...
    function cost = smoothL1(p)
        % now trying Huber loss
        delta = 0.05;
        cx = p(1); cy = p(2); ly = p(3); ry = p(4);
        cost = zeros(size(x));
        for ii = 1:length(x)
            cost(ii) = y(ii) - evalPWL(x(ii),x(1),ly,cx,cy,x(end),ry);
        end
        idx1 = find(abs(cost)>delta);
        idx2 = find(abs(cost)<=delta);
        cost(idx1) = delta*(abs(cost(idx1))-delta/2);
        cost(idx2) = cost(idx2).^2/2;
%         % for debugging, to see how many still used squares
        length(idx2) 
        length(cost)
        cost = sum(cost);
    end
end

