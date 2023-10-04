function frames = CAadjust(folder, field, times, fieldNames)
%creates structure of frames of a particular field after adjusting for movement 
%usage: CAadjust('./190320/', 33, 1:120, {'bf', 'gfp', 'cherry'});

j = field;

%intialize struct
frames = struct('field',{},fieldNames{2},{},fieldNames{3},{});
Ofield = struct('field', {}); 

%get image names and xy/t/c parameters
[xy, t, c] = imdata(dir(strcat(folder,'*xy*')));
imag = dir(strcat(folder,'*xy*'));

%Set up geometric transformation (to correct for movement between time
%points)
[opt, ~] = imregconfig('monomodal');
met = registration.metric.MattesMutualInformation;
opt.MaximumIterations = 300;

%loop on time
for i = times
    
    i

    %get images
    field = imread(strcat(folder,imag(xy==j&t==i&c==1).name));
    if max(c)>1
        field2 = imread(strcat(folder,imag(xy==j&t==i&c==2).name));
    end
    if max(c)>2
        field3 = imread(strcat(folder,imag(xy==j&t==i&c==3).name));
    end
    
    Ofield(i).field = field;
    
    %from 2cd time point start to transform 
    if (i>1)
        if (i>2)
            %from 3rd timepoint need to accumulate transformations
            oldtform = tform;
            tform = imregtform(field, Ofield(i-1).field, 'translation',opt,met,'DisplayOptimization', 0);
            temp = oldtform.T + tform.T;
            tform = affine2d([1 0 0; 0 1 0; temp(3) temp(6) 1]);
        else
            %2cd time point only transform
            tform = imregtform(field, Ofield(i-1).field,'translation',opt,met,'DisplayOptimization', 0);
        end
        %transform the field
        field = imwarp(field,tform,'OutputView',imref2d(size(Ofield(1).field)));
        if max(c)>1
            field2 = imwarp(field2,tform,'OutputView',imref2d(size(Ofield(1).field)));
        end
        if max(c)>2
            field3 = imwarp(field3,tform,'OutputView',imref2d(size(Ofield(1).field)));
        end
    end
            
        
    %save fields
    frames(i).field = field;
    if max(c)>1
        frames(i).(fieldNames{2}) = field2;
    end
    if max(c)>2
        frames(i).(fieldNames{3}) = field3;
    end
    
    %clear some variables
    clear field field2 field3
    
end