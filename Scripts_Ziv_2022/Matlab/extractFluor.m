function [fluorStruct] = extractFluor(fluorimage,rg)
%extractFluor extract pixels and properties of image based on a set of
%pixelIdxList 
%   Detailed explanation goes here
fluorStruct = struct('values',{},'mean', {} ,'std', {},'max', {} ,'top100', {} );
%tic
for i = 1:size(rg,1)
    temp = fluorimage(rg(i).PixelIdxList);
    fluorStruct(i).values = temp;
    fluorStruct(i).mean = mean2(temp);
    fluorStruct(i).std = std2(temp);
    fluorStruct(i).max = max(max(temp));
    fluorStruct(i).top300 = mean(maxk(temp,300));
end
%toc

end

