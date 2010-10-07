%%%%%%%%%
%   This program continues with previous steps. It allows user to pick
%   two points on ruler sticker in the image and convert length in pixels to
%   values in centimeters.
%                                                                        
%   ---Xiang Mao modified on Apr 4th, 2010---
%%%%%%%

clear all
close all
clc

ROOTPATH = 'C:\Documents and Settings\Xiang Mao\My Documents\MATLAB\temporary save 04-Apr-2010\'; % the folder where subfolder for each image been saved
subgroup = 's*.jpg';
sFolders = dir(fullfile([ROOTPATH],subgroup));
nff = size(sFolders,1)

    newPATH = ['C:\Documents and Settings\Xiang Mao\My Documents\MATLAB\temp_ruler_' date '\']; % the folder for saving data
    imagePATH = ['C:\Documents and Settings\Xiang Mao\My Documents\MATLAB\']; % the path for all original images
    mkdir(newPATH);

    rulertxt = ['ruler_' date '.txt'];  % the text file for save data
    fid_w = fopen([newPATH rulertxt],'w');
    fprintf(fid_w,'%s\t\n','ImageName/Pixel Per Centimeter/x(1)/x(2)/y(1)/y(2)/Orignal Ruler unit');
    
for ii = 1:nff
    ii
    iName = sFolders(ii).name
    iName_s = iName(1:(find(iName(:)=='.')-1));

    [Y]= imread([imagePATH iName]);
    imshow(Y); impixelinfo;
    
    disp('zoom in and then hit "enter"')
    pause
    disp('click on two points on the ruler to define one centimeter (or one inch), and then hit "enter"\n')
    [xr,yr]=ginput(2);
    zoom out
    xd = abs(xr(1)-xr(2));
    yd = abs(yr(1)-yr(2));
    cd = (xd^2+yd^2)^.5;
    
    ruler = '';

    while isempty(ruler)
        ut = input('Is the unit in inches? (enter "inch" if yes, otherwise enter "cm")\n','s');
        if strcmp(ut,'inch')
        ruler = 'inches'
        ppc = cd/2.54
        elseif strcmp(ut,'cm')
        ruler = 'cm'
        ppc = cd
        end
    end
    
    fprintf(fid_w,'%s\t%g\t%g\t%g\t%g\t%g\t%s\t\n',iName, ppc, xr(1),xr(2),yr(1),yr(2),ruler);
    close all
end

fclose(fid_w);
