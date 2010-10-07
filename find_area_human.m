clear all % calculation of wound area using boundary file...
close all

basepath='J:\DATA\McGuire\Traced';
filepath='\';

files = dir([basepath filepath]);
filenames = char(files.name);
filenames = filenames(3:size(filenames,1),:);

q = 1;
while q <= size(filenames,1)
    filename1 = filenames(q ,:);
    filename2 = filename1(filename1 ~= ' ');
    filename=filename2; % Use a file other than the most recent file in 'filepath'

    [Y,Ymap]= imread([basepath filepath filename]);


    imshow(Y),%('Original Image');
    impixelinfo;

    title('Select two points to reduce size of image (include the wound and some of the ruler)','FontSize',14)
    [x,y] = ginput(2); % pick two points from the image

    Ysmall=Y(y(1):y(2),x(1):x(2),:); %Y4 is the selected small picture.
    Y4=Ysmall(:,:,2); %use blue channel
    % Y4=rgb2gray(Ysmall);
    close all
    imshow(Ysmall)
%     figure
%     imshow(Y4)
%     title('Select a Wound Reference Region','FontSize',14)
%     [x,y] = ginput(2); % pick two points from the image
%     x=round(x);
%     y=round(y);
%     woundregion = Y4(y(1):y(2),x(1):x(2),:);
%     meanwound = mean(mean(woundregion));
% 
%     title('Select a Non-Wound Reference Region','FontSize',14)
%     [x,y] = ginput(2); % pick two points from the image
%     x=round(x);
%     y=round(y);
%     nonwoundregion = Y4(y(1):y(2),x(1):x(2),:);
%     meannonwound = mean(mean(nonwoundregion));

    figure
%     thresh=mean([meannonwound meanwound])/255;
    thresh=245/255;
    BW = im2bw(Y4,thresh); %Create binary image
    BW=1-BW;
    imshow(BW);
    impixelinfo;

    title('Zoom in on small region of boundary and then press enter','FontSize',14)
    zoom on
    pause
    title('Click a white point on the boundary','FontSize',14)
    [col,row] = ginput(1);
    col=round(col);
    row=round(row);
    boundary = bwtraceboundary(BW,[row, col],'E',8, 50000,'clockwise');

    imshow(Ysmall)
    hold on;
    plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);

    [mm, nn]=max(boundary);
    [mn, ni]=min(boundary);
    s=0;
    for i=mn(1):1:mm(1); 
        %ma=max(boundary(find(boundary==i),2)); 
        %mi=min(boundary(find(boundary==i),2));
        c=(find(boundary(:,1)==i));

        ma=max(boundary(c,2));
        mi=min(boundary(c,2));
        s=s + ma-mi;
        %hold on
    end
    s   % show the area of this wound in pixels
    
    button = questdlg('Is the boundary accurate?','Boundary Check','Yes','No, try again','No, skip to next picture','Yes');
    if size(button,2) == size('No, try again',2)
        close all
        continue
    elseif size(button,2) == size('No, skip to next picture',2)
        'skip'
        q = q+1;
        close all
        continue;
    end

    % figure % creats a figure window, so that previous figure window is kept
    imshow(Y(y(1):y(2),x(1):x(2),:))
    hold on
    plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);
    %plot(boundary(:,2)+x(1),boundary(:,1)+y(1),'g','LineWidth',2);

    title('Zoom in on 1cm region of ruler and then press enter','FontSize',14)
    zoom out
    zoom on
    pause

    title('Click 1 point of ruler 1=','FontSize',14)
    [x1,y1] = ginput(1);
    title('Click 1 point of ruler 2=','FontSize',14)
    [x2,y2] = ginput(1);

    xd = round(abs(x1-x2));
    yd = round(abs(y1-y2));
    cd = round((xd^2+yd^2)^.5);
    r = 10/cd;
    pixelpermm = cd/10
    areamm2 = s*r^2
    title(['Wound Area = ' num2str(areamm2)],'FontSize',14)
    zoom out

    fid = fopen([basepath '181421.txt'],'a');
    fprintf(fid,'%g\t%g\t%g\t%g\t%s\t\n',areamm2,s,pixelpermm,q,filename);
    fclose(fid)

    q = q+1;
    close all
end



