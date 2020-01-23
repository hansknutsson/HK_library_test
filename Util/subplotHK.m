function h = subplotHK(m,n,p,d,s);

%   TIGHTSUBPLOT Create axis in tiled positions tighter than SUBPLOT
%
%   H = SUBPLOTHK(M,N,P) works just like SUBPLOT but with axis tighter together.
%
%   H = SUBPLOTHK(M,N,P,D,S) sets the space between the axes to D (default 0.16)
%   relative to the axis size and the global full plot scale to S (defautlt 0.91).
%
%   For example
%
%   SUBPLOTHK(M,N,P,0.2,0.9) 
%       
%   sets the space to 20% of the axix size, and the global plot scale to 0.9. 
%
%   See also SUBPLOT and PLOT
%
%   Hans Knutsson 2017 (edited Magnus Borga's 'tightsubplot' 1999)

if nargin < 4;  d = 0.16;  end
if nargin < 5;  s = 0.91;  end

% Calculate width and height for the subplots
width   = 1/(n + (n-1)*d);
height  = 1/(m + (m-1)*d);

% Calculate position for the current subplot
row    = floor((p-1)/n) + 1;
left   = mod(p-1,n)*width*(1+d);
bottom = 1 - height*(row + (row-1)*d);
%bottom = 1 - (row*height + (row-1)*d*height)

% shrink towards center
left   = 0.5 + s*(left   - 0.5);
bottom = 0.5 + s*(bottom - 0.5);
width  = s*width;
height = s*height;

% Call subplot
ax = subplot('position',[left bottom width height]);

% return identifier, if requested:
if(nargout > 0)
    h = ax;
end
