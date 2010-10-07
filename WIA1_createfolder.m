%%%%%%%%%%%%%%
%   This program opens images in a folder one by one, creates a folder
%   for each image and allows user to choose ROI and polygon around wound
%   boundary. The croped ROI image and polygon position, as well
%   other infomation, will be saved in the corresponding folder.
%
%   ---Xiang Mao, modified on Apr 4th, 2010---
%%%%%%%%%%%%%%

%%% get file
clear all
close all
clc

rootpath = 'C:\Documents and Settings\Xiang Mao\My Documents\MATLAB\'; %the path where images to be analyzed located
ImF = dir(fullfile([rootpath],'s*.jpg')); % the common part from titles of these images
n_ImF = size(ImF, 1);

for rr = 1:n_ImF
    clear polygon*, clear Y*, clear ROI*;
    imagepath = rootpath;
    filename = ImF(rr).name;
    
    if not(ImF(rr).isdir);
        image_loaded = strcat(imagepath, filename);
        filename = strrep(filename,'JPG','jpg');
        filename_s  = filename(1:(find(filename(:)=='.')-1));
        savepath  = [imagepath '\temporary save ',date,'\'];
        mName = char(mfilename) % the current running m file
        plotpath = [savepath filename '\'];
        mkdir(plotpath);
        txtname = [mName,'_',date,'.txt']; % only save combined result, save to savepath
        infoldertxt = [filename_s '.txt']; % save all information, save to specified folder for this image
        
        Y = imread(image_loaded);
        figure(1),imshow(Y); impixelinfo; title(filename);
        
        %%% crop ROI, make sure wound in in center
        s_cr = '';
        while isempty(s_cr)
            disp('Step 1: crop a region of interest (ROI) for wound');
            disp('Click two points to define vertexes of a rectangluar that includes the entired wound region');
            figure(1);imshow(Y); impixelinfo; title(filename);
            [x,y] = ginput(2);
            yc(1)=min(y);
            yc(2)=max(y);
            xc(1)=min(x);
            xc(2)=max(x);
            ROI = Y(yc(1):yc(2),xc(1):xc(2),:);
            
            clear x, clear y;
            figure(1),clf,imshow(ROI);impixelinfo;title([filename '\_ROI']);
            [M,N,O] = size(ROI);
            
            %%% get center region position
            clear ct
            ct = [round((yc(2)-yc(1))/2) round((xc(2)-xc(1))/2)];
            cty = [ct(1)-round(M*0.10) ct(1)+round(M*0.10)];
            ctx = [ct(2)-round(N*0.10) ct(2)+round(N*0.10)];
            hold on
            plot(ct(2),ct(1),'b*');
            hold on,
            plot(ctx(1):ctx(1),cty(1):cty(2),'g.',ctx(2):ctx(2),cty(1):cty(2),'g.',...
                ctx(1):ctx(2),cty(1):cty(1),'g.',ctx(1):ctx(2),cty(2):cty(2),'g.');
            
            disp('If satisfied with the centerbox postion, hit any letter follow by "enter" to continue');
            s_cr = input('Otherwise, hit "enter" if you want to redo this step : \n','s');
            
            imwrite(ROI, [ plotpath filename_s '_ROI.jpg']);
            
            close all
            
        end
        
        %%% draw a polygon around wound region to create a skin model
        disp('Step 2: draw a polygon around wound');
        disp('The border of polygon should not cross any wound region.');
        disp('Finish drawing the polygon by a double click.');
        polygon = roipoly(ROI); %%% a logical mask
        
        %%% save workspace
        clear Y;
        save([plotpath filename_s]);
        close all
        
    end
    
end
