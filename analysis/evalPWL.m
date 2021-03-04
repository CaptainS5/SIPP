% this function is called by changeDetect.m
% it's part of the piecewise linear fit
% 26-04-2018    XW added some comments
function y = evalPWL(x,lx,ly,cx,cy,rx,ry)
warning off;
% this calculates the y value of a two-line function as plotted in
% changeDetect.m line 4-7
if x <= cx
    alpha = (x - cx)/(lx - cx);
    y = alpha*ly + (1-alpha)*cy; % it's the simplified form of cy-alpha*(cy-ly)
else
    alpha = (x - cx)/(rx - cx);
    y = alpha*ry + (1-alpha)*cy; % it's the simplified form of cy+alpha*(ry-cy)
end
end
    
