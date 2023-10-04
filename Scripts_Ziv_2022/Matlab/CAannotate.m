
function [cords, cells] = CAannotate(frames, cords, cells, cell, movie)
%allows you to annotate cells 
%usage: CAannotate(frames,[],{},'1.1','~/Desktop/CellAsic3/190320.4/anno190320.4.avi');

%intialise cord
cord = zeros(size(frames,2),1,2);

%starting time
sta = input('What time do you want to start?');
%ending time
en = input('What time do you want to end?');

%annotate
for i = sta:en
    %create frame
    %frame = im2uint8(cat(3,frames(i).field+imadjust(frames(i).mcherry,[250/2^16,5000/2^16],[0,1]),frames(i).field+imadjust(frames(i).gfp,[700/2^16,2000/2^16],[0,1]),frames(i).field));
       
    frame = im2uint8(cat(3,frames(i).field,frames(i).field+imadjust(frames(i).gfp,[700/2^16,4000/2^16],[0,1]),frames(i).field));
    %if there are already cells, show them while annotating 
    if size(cords,2)>0
        im = insertText(frame,[cords(i,:,1);cords(i,:,2)].',cells,'BoxOpacity',0,'AnchorPoint','Center');
        %im = insertText(frame,[cords(i,:,1);cords(i,:,2)].',strsplit(num2str(1:size(cords,2))),'BoxOpacity',0,'AnchorPoint','Center');
        imshow(im)
    else
        imshow(frame)
    end
    i
    [x, y] = ginput(1);
    if ~isempty(x)
        cord(i,1,1) = x;
        cord(i,1,2) = y;
    end
end

%Add cell
cords = cat(2,cords,cord);
cells = cat(1, cells, cell);

mov = VideoWriter(movie,'MPEG-4');
mov.FrameRate = 3;
open(mov)
for i = 1:size(frames,2)
    %frame = im2uint8(cat(3,frames(i).field+imadjust(frames(i).mcherry,[250/2^16,5000/2^16],[0,1]),frames(i).field+imadjust(frames(i).gfp,[700/2^16,2000/2^16],[0,1]),frames(i).field));
     
    frame = im2uint8(cat(3,frames(i).field,frames(i).field+imadjust(frames(i).gfp,[700/2^16,4000/2^16],[0,1]),frames(i).field));
    im = insertText(frame,[cords(i,:,1);cords(i,:,2)].',cells,'BoxOpacity',0,'AnchorPoint','Center');
    writeVideo(mov,im);
end
close(mov)
    
    







