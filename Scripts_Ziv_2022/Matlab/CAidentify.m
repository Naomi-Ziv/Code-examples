
function rprop = CAidentify(frames)
%identifies cells in frames of a particular field 
%usage: CAidentify(frames);

%intialize struct
names = fieldnames(frames);
name2 = names{2};
name3 = names{3};
rprop = struct('bwim', {}, 'rg', {}, name2, {}, name3, {});

%which colors empty
c = 1;
if ~isempty(frames(1).(names{2}))
    c = 2;
end
if ~isempty(frames(1).(names{3}))
    c = 3;
end


%loop on frames (times)
parfor i = 1:size(frames,2)
    i
    tempimage = frames(i).field;
    
    %bottom hat filtering - emphsize dark rings around cells
    botimage = imbothat(tempimage,strel('disk',5));
    %binarize
    bwimage = imbinarize(botimage,0.05);
    bwim = bwimage;
  
    %top hat filtering - emphsize light pixels
    %trying to break up or seperate stuff inside cell 
    topimage = imtophat(tempimage,strel('disk',5));
    tpimage = imbinarize(topimage,0.08);
    tpimage = bwmorph(tpimage,'skel',inf);
    tpimage = bwmorph(tpimage,'dilate');
    bwim(tpimage) = 0;
    
    %get rid of small islands inside cells (depends on objective)
    regprop = regionprops(bwconncomp(bwim,4),'Area','PixelIdxList','Perimeter','Eccentricity');
    bwim(cat(1,regprop([regprop.Perimeter]./[regprop.Area] > .9).PixelIdxList)) = 0;
    bwim(cat(1,regprop([regprop.Eccentricity] < .9 & [regprop.Area] < 500).PixelIdxList)) = 0;
    %get rid of small islands inside cells (depends on objective)
    bwim = bwmorph(bwim,'close');
    regprop = regionprops(bwim,'Area','PixelIdxList');
    bwim(cat(1,regprop([regprop.Area] < 100).PixelIdxList)) = 0;
    
    for int=1:3
        %trying to close up the cells
        %dilate
        bwim = bwmorph(bwim,'dilate');
        bwim = bwmorph(bwim,'thicken',4);
        %try to close
        bwim = bwmorph(bwim,'close', inf);
        %thin outlines of cells 
        bwim = bwmorph(bwim,'skel',inf);
    end
    
    %try to get cells to close
    bwim = bwmorph(bwim,'spur', inf);
    bwim = bwmorph(bwim,'close', inf);
    bwim = bwmorph(bwim,'diag');
    %invert
    bwim = bwmorph(~bwim,'thin');
    %tempimage(~bwim) = Inf;
    %imshow(tempimage)
    
    
    %save b/w image
    rprop(i).bwim = bwim;
    
    %get info
    rg = regionprops(bwim,'Area','PixelIdxList','Centroid','Perimeter','BoundingBox');
    rprop(i).rg = rg;
    
    %get fluorescent properties
    if c == 2
        rprop(i).(name2) = extractFluor(frames(i).(name2), rg);
    end
    if c == 3
        rprop(i).(name2) = extractFluor(frames(i).(name2), rg);
        rprop(i).(name3) = extractFluor(frames(i).(name3), rg);
    end
 
end