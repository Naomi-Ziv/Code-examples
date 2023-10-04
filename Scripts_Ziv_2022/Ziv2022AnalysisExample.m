%Analysis for paper: Multiple molecular events underlie stochastic switching between two heritable cell states in a eukaryotic system
%Plos Biology, 2022
%Naomi Ziv

%3 main experiments used 190320 190404 201019
%The following code was copied from a "working script" that went through
%various analysis steps, loop iterating numbers may not represent all
%fields that were analysed. Certain lines were repeatedly run with slight changes (for example: cell name)  

%also not included is various code for creating movies to see/check data/analysis visually

%General workflow included:
%1) Using CAadjust to line up frames - this creates the "frames" files which contain all the data from images 
%2) Using CAidentify for identifing cells - this creates the "rprop" files
%which contain all the cell properties data
%3) Using CAannotate recursively to annotate cells over time - this creates "cords" files
%4) using CAtrace to extract data on tracked cells over time - creates the
%various files in the "files" folder


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%190320

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Adjusting images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /Volumes/Leeuwenhoek/

exp = '190320';
home = '/Volumes/Newton/190320.2/frames';

parfor t = 1:3  % j = 72:77 %j = 31:71
    
    temp = [5 10 82];
    j = temp(t);
    j
    
    name = strcat(exp,'.',num2str(j));
    
    frames = CAadjust(strcat('./',exp,'/'), j, 1:120, {'bf', 'gfp', 'mcherry'});
    
    parsave(frames,home,name,'frames.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Identifing cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vec = [31 33 35 37 40 43 45 48 50 54 61 64 66 67 69 74];
%vec = [61 66 67 69 74];
%vec = [33 37 40 45 64 66 67 69];
%vec = [36 38 39 77];
vec = [80];

parfor i = 1:size(vec,2)
    j = vec(i);
    j
    
    framename = strcat(home,'/frames.190320.',num2str(j),'.mat');  
    frames = load(framename);
    frames = frames.frames;
    rprop = CAidentify(frames);
    
    name = strcat(exp,'.',num2str(j));
    parsave2(rprop, '/Volumes/Newton/190320.2/rprop',name, 'rprop.')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Annotating cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '190320';
home = '/Volumes/Newton/190320.2/';

j=33
name = strcat(exp,'.',num2str(j));
load(strcat(home,'frames/frames.',name,'.mat'))

[cords, cells] = CAannotate(frames,[],{},'0',strcat('/Volumes/Newton/190320.2/movies2/cords.',name,'.mp4'));
[cords, cells] = CAannotate(frames,cords,cells,'1',strcat('/Volumes/Newton/190320.2/movies2/cords.',name,'.mp4'));

save(strcat(home,'/cords/cords.',name,'.mat'),'cords','cells','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Extracting data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '190320';
home = '/Volumes/Newton/190320.2/'

for j = [80]
   
    name = strcat(exp,'.',num2str(j));
    
    cordname = strcat(home,'cords/cords.190320.',num2str(j),'.mat');
    rpropname = strcat(home,'rprop/rprop.190320.',num2str(j),'.mat');
    
    rprop = load(rpropname);
    rprop = rprop.rprop;
    cords = load(cordname);
    cells = cords.cells;
    cords = cords.cords;

    traces = CAtrace(rprop,cords,'gfp');
    
    temp=traces.mean;save(strcat(home,'/files/','Means.',name,'.txt'),'temp','-ascii')
    temp=traces.area;save(strcat(home,'/files/','Areas.',name,'.txt'),'temp','-ascii')
    temp=traces.max;save(strcat(home,'/files/','Maxs.',name,'.txt'),'temp','-ascii')
    temp=traces.top300;save(strcat(home,'/files/','Tops.',name,'.txt'),'temp','-ascii')

    writetable(table(cells),strcat(home,'/files/','Cells.',name,'.txt'),'WriteVariableNames',false)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%190404

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Adjusting images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /Volumes/Leeuwenhoek/

exp = '190404';
home = '/Volumes/Newton/190404.2/frames';
vec = [96 106 99 123 126 144];
parfor i = 1:size(vec,2)
    j = vec(i);
    j
    
    name = strcat(exp,'.',num2str(j));
    
    frames = CAadjust(strcat('./',exp,'/'), j, 1:120, {'bf', 'gfp', 'mcherry'});
    
    parsave(frames,home,name,'frames.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Identifing cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '190404';
home = '/Volumes/Newton/190404.2/frames';
vec = [114 155 96 106 99 123 126 144];
parfor i = 1:size(vec,2)
    j = vec(i);
    j
    
    framename = strcat(home,'/frames.190404.',num2str(j),'.mat');  
    frames = load(framename);
    frames = frames.frames;
    rprop = CAidentify(frames);
    
    name = strcat(exp,'.',num2str(j));
    parsave2(rprop, '/Volumes/Newton/190404.2/rprop',name, 'rprop.')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Annotating cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '190404';
home = '/Volumes/Newton/190404.2/';

j=106
name = strcat(exp,'.',num2str(j));
load(strcat(home,'frames/frames.',name,'.mat'))

[cords, cells] = CAannotate(frames,[],{},'0',strcat('/Volumes/Newton/190404.2/movies2/cords.',name,'.mp4'));
[cords, cells] = CAannotate(frames,cords,cells,'1',strcat('/Volumes/Newton/190404.2/movies2/cords.',name,'.mp4'));

save(strcat(home,'/cords/cords.',name,'.mat'),'cords','cells','-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Extracting data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '190404';
home = '/Volumes/Newton/190404.2/';

for j = [106]
   
    name = strcat(exp,'.',num2str(j));
    
    cordname = strcat(home,'cords/cords.190404.',num2str(j),'.mat');
    rpropname = strcat(home,'rprop/rprop.190404.',num2str(j),'.mat');
    
    rprop = load(rpropname);
    rprop = rprop.rprop;
    cords = load(cordname);
    cells = cords.cells;
    cords = cords.cords;

    traces = CAtrace(rprop,cords,'gfp');
    
    temp=traces.mean;save(strcat(home,'/files/','Means.',name,'.txt'),'temp','-ascii')
    temp=traces.area;save(strcat(home,'/files/','Areas.',name,'.txt'),'temp','-ascii')
    temp=traces.max;save(strcat(home,'/files/','Maxs.',name,'.txt'),'temp','-ascii')
    temp=traces.top300;save(strcat(home,'/files/','Tops.',name,'.txt'),'temp','-ascii')

    writetable(table(cells),strcat(home,'/files/','Cells.',name,'.txt'),'WriteVariableNames',false)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%201019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Adjusting images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /Volumes/Newton/

exp = '201019';
home = '/Volumes/Newton/201019.2/frames';

%vec = [3   4   5  11  15  27  30  34  41  47  54  62  63  67  74  83  84  98 101 112 149 155 156 168 169 170 172 180 184 193 197 201 211 213];
%vec = [132];

vec = [123 127 142];
parfor i = 1:size(vec,2)
    
    j = vec(i);
    j
    
    name = strcat(exp,'.',num2str(j));
    
    frames = CAadjust(strcat('./',exp,'/'), j, 1:120, {'bf', 'gfp', 'mcherry'});
    
    parsave(frames,home,name,'frames.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Identifing cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /Volumes/Newton/
exp = '201019';
home = '/Volumes/Newton/201019.2/frames';
%vec = [5 41 156 184];
vec = [123 127 142];
parfor i = 1:size(vec,2)
    j = vec(i);
    j
    
    framename = strcat(home,'/frames.201019.',num2str(j),'.mat');  
    frames = load(framename);
    frames = frames.frames;
    rprop = CAidentify(frames);
    
    name = strcat(exp,'.',num2str(j));
    parsave2(rprop, '/Volumes/Newton/201019.2/rprop',name, 'rprop.')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Annotating cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '201019';
home = '/Volumes/Newton/201019.2/';

j=127
name = strcat(exp,'.',num2str(j));
load(strcat(home,'frames/frames.',name,'.mat'))

[cords, cells] = CAannotate(frames,[],{},'0',strcat('/Volumes/Newton/201019.2/movies2/cords.',name,'.mp4'));
[cords, cells] = CAannotate(frames,cords,cells,'1',strcat('/Volumes/Newton/201019.2/movies2/cords.',name,'.mp4'));

save(strcat(home,'/cords/cords.',name,'.mat'),'cords','cells','-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Extracting data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp = '201019';
home = '/Volumes/Newton/201019.2/'

for j = [127]
   
    name = strcat(exp,'.',num2str(j));
    
    cordname = strcat(home,'cords/cords.201019.',num2str(j),'.mat');
    rpropname = strcat(home,'rprop/rprop.201019.',num2str(j),'.mat');
    
    rprop = load(rpropname);
    rprop = rprop.rprop;
    cords = load(cordname);
    cells = cords.cells;
    cords = cords.cords;

    traces = CAtrace(rprop,cords,'gfp');
    
    temp=traces.mean;save(strcat(home,'/files/','Means.',name,'.txt'),'temp','-ascii')
    temp=traces.area;save(strcat(home,'/files/','Areas.',name,'.txt'),'temp','-ascii')
    temp=traces.max;save(strcat(home,'/files/','Maxs.',name,'.txt'),'temp','-ascii')
    temp=traces.top300;save(strcat(home,'/files/','Tops.',name,'.txt'),'temp','-ascii')

    writetable(table(cells),strcat(home,'/files/','Cells.',name,'.txt'),'WriteVariableNames',false)

end




