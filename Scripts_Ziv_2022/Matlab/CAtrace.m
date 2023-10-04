function [traces] = CAtrace(rprop, cords, fluor)
%extracts information for cells over time based on set of coordinates
%usage: CAtrace(frames,cords, 'gfp');

%%fix to be NaN to begin with
cords(cords==0) = NaN;
%intialize
z = zeros(size(cords,1),size(cords,2));
traces = struct('area',z,'mean', z, 'max', z ,'top300', z );

%get cell identifiers
for i = 1:size(rprop,2)
    
    %get linear index
    label = bwlabel(rprop(i).bwim);
    idx = sub2ind(size(label),(cords(i,:,2)),(cords(i,:,1)));
    idx2 = zeros(size(idx));
    idx2(~isnan(idx)) = label(int64(idx(~isnan(idx))));
    idx3 = idx2(idx2~=0);
    
    %get info
    if ~isempty(idx3)
        traces.area(i,idx2~=0) = [rprop(i).rg(idx3).Area]; 
        traces.mean(i,idx2~=0) = [rprop(i).(fluor)(idx3).mean];
        traces.max(i,idx2~=0) = [rprop(i).(fluor)(idx3).max];
        traces.top300(i,idx2~=0) = [rprop(i).(fluor)(idx3).top300];
    end
    
end

