function [xy, t, c] = imdata(directory)
% returns xy, t (time), and c (color) data about a series of tiff files

d = directory;
xy = 1:length(d);
t = 1:length(d);
col = 1:length(d);

n = d(1).name;
a = 1:4;
if ~isempty(regexp(n , 'xy0', 'once'))
    a(1) = regexp(n , 'xy0');
elseif ~isempty(regexp(n , 'xy1', 'once'))
	a(1) = regexp(n , 'xy1');
elseif ~isempty(regexp(n , 'xy2', 'once'))
  a(1) = regexp(n , 'xy2');
elseif ~isempty(regexp(n , 'xy3', 'once'))
  a(1) = regexp(n , 'xy3');
elseif ~isempty(regexp(n , 'xy4', 'once'))
  a(1) = regexp(n , 'xy4');
elseif ~isempty(regexp(n , 'xy5', 'once'))
  a(1) = regexp(n , 'xy5');
elseif ~isempty(regexp(n , 'xy6', 'once'))
  a(1) = regexp(n , 'xy6');
elseif ~isempty(regexp(n , 'xy7', 'once'))
  a(1) = regexp(n , 'xy7');
elseif ~isempty(regexp(n , 'xy8', 'once'))
  a(1) = regexp(n , 'xy8');
elseif ~isempty(regexp(n , 'xy9', 'once'))
  a(1) = regexp(n , 'xy9');
else
    a(1) = 0;
end
if ~isempty(regexp(n , 't0', 'once'))
    a(2) = regexp(n , 't0');    
elseif ~isempty(regexp(n , 't1', 'once'))
	a(2) = regexp(n , 't1'); 
else
    a(2) = 0;
end        
if ~isempty(regexp(n , 'c1', 'once'))
    a(3) = regexp(n , 'c1');    
else
    a(3) = 0;
end        
if ~isempty(regexp(n , '.tif', 'once'))
    a(4) = regexp(n , '.tif');
elseif    ~isempty(regexp(n , '.jpg', 'once'))
    a(4) = regexp(n , '.jpg');
else
    a(4) = 0;
end

if a(1)
    x = a - a(1);
    x(find(x<1)) = 1000;
    [C,I] = min(x);
    cxy = a(1)+2:a(I)-1;
else
    cxy = [];
    xy = [];
end

if a(2)
    x = a - a(2);
    x(find(x<1)) = 1000;
    [C,I] = min(x);
    ct = a(2)+1:a(I)-1;
else
    ct = [];
    t = [];
end


if a(3)
    x = a - a(3);
    x(find(x<1)) = 1000;
    [C,I] = min(x);
    cc = a(3)+1:a(I)-1;
else
    cc = [];
    col = [];
end

for i = 1:length(d)
    n = d(i).name;
    if ~isempty(cxy)
        xy(i) = str2num(n(cxy));
    end
    if ~isempty(ct)
        t(i) = str2num(n(ct));
    end
    if ~isempty(cc)
        col(i) = str2num(n(cc));
    end
end
c = col;
if isempty(c)
    c = 1:length(xy);
    c(:) = 1;
end