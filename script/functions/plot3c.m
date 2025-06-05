function h = plot3c(x,y,z,c,symbol)
% PLOT3C generates a 3-D plot using the colormap for color
%
% PLOT3C(X,Y,Z) produces a 3-D plot using the data defined
% by X, Y, and Z.  Note the X, Y, and Z must be a vector. 
% This function does not support matrices as input.
%
% PLOT3C(X,Y,Z,'SYMBOL') produces a 3-D plot using the symbol
% defined by the string SYMBOL.  Valid symbols are:
%
% [ {-} | -- | : | -. | none | + | o | * | . | x | square | 
% |diamond | v | ^ | > | < | pentagram | hexagram | {none} ]  
%
%
% PLOT3C(X,Y,Z,C) maps that data in C to the colormap.  This
% is the CData property.  See the description of CData under
% the surface command in the MATLAB Reference Guide for more
% information.
%
% PLOT3C(X,Y,Z,C,'SYMBOL') lets you specify everything.
%
% H = PLOT3C(...) returns the handle to the surface plot.
%
% EXAMPLE:
%
%    t = linspace(0,2*pi,40);
%    x = sin(t);
%    y = cos(t);
%    z = t.*x.^2;
%    h = plot3c(x,y,z,t,'x');
%
% SEE ALSO: SURFACE, SURF, MESH
% DISCLAIMER:  This M-file was written as an example of how
%              to generate 3-D line plots that use the colormap.
%              It has not tested by The MathWorks, Inc., and
%              therefore, it is not supported.  Please feel free
%              to edit this file to suit your needs.

% modified by Nathan Briggs
% 11 March, 2013
% duplicates points if length(x) == 1
% 11 January, 2011
% eliminates NANs
% 17 September, 2009
% forces c to be a double precision floating point variable, allowing for
% input variables of format single or int

% Written by John L. Galenski III
% Copyright (c) 1995 by The MathWorks, Inc.
% All Rights Reserved

% Parse the inputs
if nargin < 3
  error('Requires 3 or more inputs.')
end

if nargin == 3   % X,Y,Z given
  c = z;
  symbol = '-';
end

if nargin == 4  % X,Y,Z, and C or SYMBOL given
  if isstr(c)   % SYMBOL given
    symbol = c;
    c = z;
  else          % C given
    symbol = '-';
  end
end

x = x(:); y = y(:); z = z(:); c = c(:);
notnan = true(size(c));
if ~isa(x,'datetime')
    notnan = notnan&~isnan(x);
end
if ~isa(y,'datetime')
    notnan = notnan&~isnan(y);
end
if ~isa(z,'datetime')
    notnan = notnan&~isnan(z);
end
if ~isa(c,'datetime')
    notnan = notnan&~isnan(c);
end
x = x(notnan); y = y(notnan); z = z(notnan); c = c(notnan);
h = nan;
if isempty(x)
    warning('x contains no data')
elseif isempty(y)
    warning('y contains no data')
elseif isempty(z)
    warning('z contains no data')
elseif isempty(c)
    warning('c contains no data')
else
    if length(x) == 1
        x = [x; x];
        y = [y; y];
        z = [z; z];
        c = [c; c];
    end
    c = double(c);
    
    % Generate the mesh plot
    X = [x(:),x(:)];
    Y = [y(:),y(:)];
    Z = [z(:),z(:)];
    C = [c(:),c(:)];
    h = mesh(X,Y,Z,C);
    if ismarker(symbol)
        set(h, 'Marker', symbol, 'LineStyle', 'none')
    else
        set(h,'LineStyle',symbol)
    end
end

function y = ismarker(symbol)
%ISMARKER determines if the specified symbol is a valid 
%marker or a valid linestyle propery value.

markers = {'+';'o';'*';'.';'x';'square';'diamond';
           'v';'^';'>';'<';'pentagram';'hexagram';'none'};
for i=1:length(markers) 
   val=strcmp(markers{i}, symbol);
   if val==1
      break
   end
end

if nargout>0
   y=val;
end

