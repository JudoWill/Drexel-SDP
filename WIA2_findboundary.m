%%%%%%%%%
%   This program continues with the first step. It performs analysis on
%   the ROI image, finds boundary on multiple color channels and differet
%   threshold coefficients. Pixels on boundaries found at under each situation
%   are saved by their coordinates in the ROI image.
%
%   ---Xiang Mao modified on Apr 4th, 2010---
%%%%%%%


clear all
close all
clc
tic

%Will:Need to make these varaibles and not hard-coded
ROOTpath = 'C:\Documents and Settings\Xiang Mao\My Documents\MATLAB\temporary save 04-Apr-2010\'; % where the folders for each image been located
folders = dir(fullfile([ROOTpath],'*.jpg'));
n_folders = length(folders);

%%

for rr = 1:n_folders
    %Will: should NEVER use clear like this ... it makes debugging
    %Will: IMPOSSIBLE
    clear matname, clear matpath,
    clear ROI*, clear BW*, clear v*, clear polygon*, clear Y*,
    clear center*, clear perim*,
    
    matname = strrep(folders(rr).name, '.jpg','');
    matpath = [ROOTpath folders(rr).name '\']; 
    %Will: '\' assumes Windows PC's not always true!
    %Will: "fullfile" or "pathsep" is a much safer alternative
    
    if folders(rr).isdir
        
        load([matpath matname]); % load the '.mat' file
        %Will: use fullpath because its safer!
        %Will: also should probably load things in as a struct to unpack
        %Will: them.  There's no way to see what variables this has
        %Will: introduced!
        
        mName = char(mfilename);
        txtname = [mName,'_',date,'.txt'];
        
        savepath = ROOTpath;
        plotpath = matpath;
        
        ROI_loaded = strcat(matpath, filename_s, '_ROI.jpg');
        %Will: use fullpath because its safer!
        
        fid2 = fopen(fullfile(plotpath,infoldertxt),'w'); % detailed infomation.
        %%%%% infoldertxt was saved in '.mat' file from previous step.
        
        ROI = imread([ROI_loaded]);
        
        perim_polygon = bwperim(polygon);
        polygon_inv = logical(1 - polygon);
        polygon3 = cat(3,polygon, polygon, polygon );
        
        se3 = strel('disk',3);
        se5 = strel('disk',5);
        se7 = strel('disk',7);
        
        ROI_R = ROI(:,:,1);
        ROI_G = ROI(:,:,2);
        ROI_B = ROI(:,:,3);
        
        srgb2lab = makecform('srgb2lab');
        ROI_Lab = applycform(im2double(ROI), srgb2lab); % convert to L*a*b*
        
        ROI_L = ROI_Lab(:,:,1);
        ROI_a = ROI_Lab(:,:,2);
        ROI_b = ROI_Lab(:,:,3);
        
        center_Lab = ROI_Lab(cty(1):cty(2), ctx(1):ctx(2),:);
        center_Lab_int = [mean2(center_Lab(:,:,1)) mean2(center_Lab(:,:,2)) mean2(center_Lab(:,:,3))];
        
        ROI_E = sqrt((ROI_Lab(:,:,1)-center_Lab_int(1)).^2+ (ROI_Lab(:,:,2)-center_Lab_int(2)).^2+(ROI_Lab(:,:,3)-center_Lab_int(3)).^2);
        ROI_Er = sqrt((ROI_Lab(:,:,1)-54).^2+ (ROI_Lab(:,:,2)-81).^2+(ROI_Lab(:,:,3)-70).^2);
        
        fod = 10;  % Factor Of Reduction, for Pixel-Color-Compare.
        
        ROI_reduc(:,:,:) = fod*round(ROI(:,:,:)/fod); %reduce colors.
        
        %% loop in channels/layers
        skinmodel_marker = 0; %%%if the code has been run on skin model (= pixel color compare) once, the marker will be turned to 1.
        COEF = linspace(0,1,11);
        
        allchannels = {'R','G','B','L','a','b','E','Er','reduc'};
        channelweight = [0, 2, 1, 1, 2,	0, 2, 0, 1];   %% Note: these weights is what been used in the paper. XM. 4/4/2010
        
        b = zeros(1,length(allchannels));
        a = char(allchannels);
        for ff = 1:length(allchannels)
            b(ff) = channelweight(ff);
        end
        channelnumber = length(allchannels);
        n = channelnumber;
        
        
        for cc = 1: 11
            cof = COEF(cc);
            totalweight = 0;
            
            fprintf(fid2,'%s\t%4.5g\t%4.5g\t%4.5g\t%4.5g\t%4.5g\t%4.5g\t\n','clock', clock);
            fprintf(fid2,'%s\t%s\t\n','path', imagepath);
            fprintf(fid2,'%s\t%s\t\n','m-file', mName);
            fprintf(fid2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',a(1,:),a(2,:),a(3,:),a(4,:),a(5,:),a(6,:),a(7,:),a(8,:),a(9,:));
            fprintf(fid2,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t\n',b(:));
            fprintf(fid2,'%s\t%g\t\n','factor of reduction',fod);
            fprintf(fid2,'%s\t\n',...
                'filename/channel/area/perimeter/level/center_ave/perimeter_ave/max/min/coef');
            
            BW_combine1 = double(false(M,N));
            area(1:n) = 0;
            perimlength(1:n) = 0;
            
            for jj  = 1: n % n is channelnumber
                channelname = char(allchannels(jj));
                fprintf(1,'%s\t%g\t\t\t\t%s\t%s\t\n','coefficient:',cof,'channel name:',channelname);
                clear layer, clear BW, clear v,
                close all
                
                eval(['layer = ','ROI_',channelname,';']);
                %Will: EVAL is the devil's minion, it should never be used
                %Will: like this!
                % % % %     figure(20),imshow(layer,[]);impixelinfo; title(channelname);
                
                if not(strcmp('reduc', channelname))
                    %%% all channels that need coefficients
                    layer = double(layer);
                    
                    perim = immultiply(perim_polygon, layer);
                    perim_list = double(nonzeros(perim));
                    perim_ave = mean(perim_list);
                    perim_dev = std(double(perim_list));
                    
                    center = layer(cty(1):cty(2), ctx(1):ctx(2),:);
                    center_list = double(nonzeros(center));
                    center_ave = mean(center_list);
                    center_dev = std(double(center_list));
                    
                    level = center_ave - cof*(center_ave - perim_ave);
                    
                    ma = max(nonzeros(immultiply(layer, polygon)));
                    mi = min(nonzeros(immultiply(layer, polygon)));
                    layerrange = ma-mi;
                    
                    
                    normlayer = (layer-mi)./layerrange;
                    normlevel = (level-mi)./layerrange;
                    BW = im2bw(normlayer, normlevel);
                    
                    if center_ave < perim_ave
                        BW = 1-BW;
                    end
                    
                    BW = immultiply(polygon, BW);
                    BW = imopen(BW,se3);
                    BW = imclose(BW,se5);
                    BW = imfill(BW,'holes');
                    
                    B = bwboundaries(BW,8,'noholes');
                    d=cellfun('length', B);
                    [max_d,k]=max(d);
                    v = B{k(1)};
                    % % %             figure(100), imshow(ROI);impixelinfo;
                    % % %             title([channelname, '\_', num2str(level), '\_',num2str(cof) ]);
                    % % %             hold on;
                    % % %             plot(v(:,2),v(:,1),'g','LineWidth',2);
                    
                    
                    eval(['v_',channelname,'_',num2str(cof*10), '= v;']); % boundary info for each channel and coefficient
                    eval(['level_',channelname,'_',num2str(cof*10), '= level;']); % real threshold
                    eval(['perim_',channelname, '_ave','= perim_ave;']); % average intensity on polygon, one value for one channel
                    eval(['perim_',channelname, '_dev','= perim_dev;']);
                    eval(['center_',channelname, '_ave','= center_ave;']); % average intensity on center, one value for one channel
                    eval(['center_',channelname, '_dev','= center_dev;']);
                    eval(['max_',channelname, '= ma;']);
                    eval(['min_',channelname, '= mi;']);
                    %Will: EVAL is the devil's minion, it should never be used
                    %Will: like this!
                    
                else   %%%% this is skin model (= pixel color compare)
                    
                    while skinmodel_marker == 0;  %% Note: so this while-end loop will only be run once when COEF loop from 0 to 11.
                        
                        polygon_inv3 = cat(3,polygon_inv,polygon_inv,polygon_inv);
                        polygon_inv_ROI = immultiply(ROI_reduc, uint8(polygon_inv3));
                        
                        modellist=reshape(polygon_inv_ROI,[M*N 3]); %reshape boundary pixels to rows containing [R G B]
                        modelunique=unique(modellist,'rows'); %unique list of color combinations on the boundary
                        
                        [x y]=find(polygon ==1); %find all pixels within the polygon
                        
                        %%%%%%%%% loop through the polygon-enclosed pixels
                        for i=1:size(x,1)  % size(x,1) is the length of x, which is the # of pixels inside polygon.
                            xi=x(i);
                            yi=y(i);
                            for j=1:size(modelunique,1)
                                if layer(xi,yi,1)==modelunique(j,1) && layer(xi,yi,2)==modelunique(j,2) && layer(xi,yi,3)==modelunique(j,3)
                                    layer(xi,yi,1)=0;
                                    layer(xi,yi,2)=0;
                                    layer(xi,yi,3)=0;
                                end  %% to eliminate pixels have the same color as skin model
                            end
                        end
                        polygon3 = cat(3,polygon,polygon,polygon); % replicate the polygon for 3 channels
                        layer(find(polygon3==0))=0; % set pixels outside the polygon to zero
                        
                        % % %             figure(30),imshow(layer);impixelinfo; title('reduc');
                        
                        BW = ones(M,N);
                        for i = 1:M
                            for j = 1:N
                                if layer(i,j,1) == 0 && layer(i,j,2) == 0 && layer(i,j,3) == 0
                                    BW(i,j) = 0;
                                else
                                    BW(i,j) = 1;
                                end
                            end
                        end   %% make a binary image,
                        
                        BW = imclose(BW, se5);
                        BW = imfill(BW);
                        
                        B = bwboundaries(BW,8,'noholes');
                        d=cellfun('length', B);
                        [max_d,k]=max(d);
                        v=B{k(1)};
                        % % %             figure(100), imshow(ROI);impixelinfo;title(channelname);
                        % % %             hold on;
                        % % %             plot(v(:,2),v(:,1),'g','LineWidth',2);
                        % % %             h = gcf;
                        % % %             savechannel = ['RGB_' channelname];
                        % % %             saveas(h,[plotpath filename_s '_' savechannel '.jpg' ]);
                        
                        %%%%skin model
                        eval(['v_',channelname, '= v;']);
                        eval(['level_',channelname, '= ','NaN;']);
                        
                        eval(['perim_',channelname, '_ave','= ','NaN;']);
                        eval(['perim_',channelname, '_dev','= ','NaN;']);
                        eval(['center_',channelname, '_ave','= ','NaN;']);
                        eval(['center_',channelname, '_dev','= ','NaN;']);
                        eval(['max_',channelname, '= ','NaN;']);
                        eval(['min_',channelname, '= ','NaN;']);
                        %Will: EVAL is the devil's minion, it should never be used
                        %Will: like this!
                        
                        skinmodel_marker = skinmodel_marker + 1;
                        
                    end %end of while
                    %%%% when skinmodel_marker >0, use previous v_reduc for v;
                    
                    %%%%skin model (= pixel color compare)
                end
                
                if strcmp(channelname,'reduc');
                    v = eval(['v_',channelname]);
                    %Will: EVAL is the devil's minion, it should never be used
                    %Will: like this!
                end
                
                lv = size(v,1); % perimeter.
                Z = false(M,N);
                
                if v==0
                    disp('v = 0');
                    disp([channelname ' layer will not be used']); % this is redundent after remove ttest
                else
                    for kk = 1:lv;
                        Z(v(kk,1),v(kk,2))=1;
                    end
                    Z = imfill(Z,'hole');
                    BW_combine1 = BW_combine1 + channelweight(jj)*Z;
                    totalweight = channelweight(jj) + totalweight;
                end
                
                area(jj) = sum(sum(Z)); % area for single channel
                perimlength(jj) = lv;
                %Will: These can easily be pre-allocated!
                
                if not(strcmp('reduc', channelname))
                    fprintf(fid2,'%s\t%s\t%g\t%g\t%8.3g\t%8.3g\t%8.3g\t%8.3g\t%8.3g\t%8.3g\t\n',...
                        filename,channelname,area(jj),perimlength(jj),level,center_ave, perim_ave,ma,mi,cof);
                else
                    fprintf(fid2,'%s\t%s\t%g\t%g\t\n',filename,channelname,area(jj),perimlength(jj));
                end
                
            end
            
            close all
            clear B, clear v
            
            %%% combine all BW images that are avaiable
            BW_combine1 = BW_combine1./totalweight;
            % % % figure(200),imshow(BW_combine1,[]);impixelinfo;title(['BW\_combine1\_', num2str(cof)] );
            % % % h = gcf;
            % % % saveas(h,[plotpath filename_s '_combine1_' num2str(cof) '.jpg']);
            
            BW_combine2 = round(BW_combine1);
            BW_combine2 = imopen(BW_combine2,se3);
            BW_combine2 = imclose(BW_combine2,se7);
            BW_combine2 = imfill(BW_combine2);
            % % % figure(201),imshow(BW_combine2,[]);impixelinfo;title('BW\_combine2');
            
            
            B = bwboundaries(BW_combine2,8,'noholes');
            d=cellfun('length', B);
            [max_d,k]=max(d);
            v = B{k(1)};
            figure(300), imshow(ROI);impixelinfo;
            hold on;
            FN = strrep(filename_s,'_', '\_');
            plot(v(:,2),v(:,1),'g','LineWidth',2);title([FN  '\_combine\_' num2str(cof) ]);
            h = gcf;
            saveas(h,[plotpath filename_s '_combine_' num2str(cof) '.jpg']);
            
            
            lv = size(v,1); % perimeter.
            Z = false(M,N);
            for kk = 1:lv;
                Z(v(kk,1),v(kk,2))=1;
            end
            eval(['v_combine_',num2str(cof*10), '= v;']);
            Z = imfill(Z,'hole');
            BW_combine = Z ;
            area_c = sum(sum(Z));
            perimlength_c = lv;
            
            fprintf(fid2,'%s\t%s\t%g\t%g\t\t\t\t\t\t\t\t\t\t\t%8.3g\t\n\n',filename,'combo',area_c,perimlength_c,cof);
            
            close all
            
        end
        
        
        fclose(fid2);
        filename
        
        close all
        %% save workspace
        
        clear ROI_*;
        clear layer;
        clear se*;
        clear v;
        clear srgb2lab;
        clear polygon_inv*;
        clear polygon3;
        clear perim_polygon;
        clear a, clear b;
        clear Z, clear BW,clear BW_combine2;
        %Will: Shouldn't clear variables like this it makes debugging
        %Will: difficult!
        
        save([plotpath filename_s '_allv']);
        %Will: Should use fullfile!
        
    end
    
end %%%end the muti-folder loop of rr
time = toc


