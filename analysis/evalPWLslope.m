% this function is called by changeDetectSlope.m
% it's part of the piecewise linear fit using slopes as parameters
% Jan-11-2021    XW modified based on evalPWL.m
function y = evalPWLslope(x,lx,slope1,cx,cy,rx,slope2)
warning off;
% this calculates the y value of a two-line function as plotted in
% changeDetectSlope.m line 5-8
ly = cy-slope1*(cx-lx);
ry = cy-slope2*(cx-rx);
if x <= cx
    alpha = (x - cx)/(lx - cx);
    y = alpha*ly + (1-alpha)*cy; % it's the simplified form of cy-alpha*(cy-ly)
else
    alpha = (x - cx)/(rx - cx);
    y = alpha*ry + (1-alpha)*cy; % it's the simplified form of cy+alpha*(ry-cy)
end
end
    
