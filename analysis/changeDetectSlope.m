% FUNCTION to find pursuit onset
% find the change in y(x), by fitting a piecewise linear model
% with one break. To restrict the value of slopes, instead of fitting
% ly and ry as in changeDetect.m, here we fit for slope1 (before cx) and slope2 (after cx)
%  ly\    
%     \  /ry 
%      \/
%       (cx,cy)
% Input:
% slopeRange should be a 2x2 matrix, column1 = slope1, column2 = slope2,
%   row1 = lower bound, row2 = upper bound.
% cxRange should have a lower bound cxRange(1) and an upper bound cxRange(2)
%
% History
% Jan-11-2021 modified from changeDetect.m by Xiuyun Wu,
% xiuyunwu5@gmail.com

function [cx,cy,slope1,slope2] = changeDetectSlope(x, y, cxRange, slopeRange)
warning off; 
% initialize parameters, p0, to reasonable values
w = ceil(length(x)/10);% small window
ly = mean(y(1:w));
ry = mean(y(end-w+1:end));
cy = mean(y(w+1:end-w));
cx = mean(x(w+1:end-w));
slope1 = (cy-ly)/(cx-x(1));
slope2 = (cy-ry)/(cx-x(end));
p0 = [cx, cy, slope1, slope2];
% minimize residual
options = optimset('Display', 'off');
p = lsqnonlin(@myfun,p0,[cxRange(1), -inf, slopeRange(1, 1), slopeRange(1, 2)],[cxRange(2), inf, slopeRange(2, 1), slopeRange(2, 2)],options);
cx = p(1); cy = p(2); slope1 = p(3); slope2 = p(4);
return;
% inner function for computing residual errors
    function residual = myfun (p)
        cx = p(1); cy = p(2); slope1 = p(3); slope2 = p(4);
        residual = zeros(size(x));
        for i = 1:length(x)
            residual(i) = y(i) - evalPWLslope(x(i),x(1),slope1,cx,cy,x(end),slope2);
        end
    end

end

