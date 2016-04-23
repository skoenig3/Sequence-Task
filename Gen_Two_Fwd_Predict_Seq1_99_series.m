% Same as code for Gen_Two_Fwd_Predict.m except different seeds and instead
% of SEQ RM-RZ it is Seq01-Seq99. Modified code on 3/4/15 by Seth
% Koenig
number_of_sequences = 2;
number_of_crosshairs = 4;
crosshair_locations = {};
overlap_index = [];
buffer = 5; %minimum distance between 2 crosshairs
maxdist = 15; %maximum distance between 2 crosshairs
buffer = buffer - 1;% since items are 1 dva grid
sizex = 25;
sizey = 19;
clear x y
[cc,rr] = meshgrid(1:sizex,1:sizey);
rand('seed',150304);
for i = 1:1000;
    valid_points = ones(sizey,sizex);
    crosshair_locations{i} = NaN(number_of_crosshairs*2,number_of_sequences);
    seq = 1;
    for cross = 1:number_of_crosshairs
        if cross == 1
            available_points = find(valid_points);
            choosen_ind = available_points(randi(length(available_points)));
            [y,x] = ind2sub([sizey,sizex],choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2)<= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            crosshair_locations{i}(2*cross-1,seq) = x;
            crosshair_locations{i}(2*cross,seq) = y;
        elseif cross == 2
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            choosen_ind = Cind(randi(length(Cind)));
            [y,x] = ind2sub(size(valid_points),choosen_ind);
            C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
            Cind = find(C);
            valid_points(Cind) = 0;
            crosshair_locations{i}(2*cross-1,seq) = x;
            crosshair_locations{i}(2*cross,seq) = y;
        else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
            %find angle of quickest scan path from previous 2 crosshairs
            dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
            dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
            angle12 = atan2d(dy12,dx12);
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            
            available_points = find(valid_points);
            dist = sqrt((rr-y).^2+(cc-x).^2);
            C = dist <= 15;
            Cind = find(C);
            Cind = intersect(Cind,available_points);
            [Cy,Cx] = ind2sub(size(valid_points),Cind);
            potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
            for pc = 1:length(Cy)
                dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                angle23 = atan2d(dy23,dx23);
                potential_angles(pc) = angle23;
            end
            potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
            angle_difference = potential_angles-angle12;
            angle_difference = abs(angle_difference);
            good_angles = find(angle_difference <=  90);
            if isempty(good_angles);
                break;
            else
                choosen_ind = good_angles(randi(length(good_angles)));
                crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                x =  Cx(choosen_ind);
                y =  Cy(choosen_ind);
                Cind = find(C);
                valid_points(Cind) = 0;
            end
        end
    end
    
    if all(~isnan(crosshair_locations{i}(:,1)))
        overlapind = rem(i,4)+1;
        overlap_index(i) = overlapind;
        seq = 2;
        switch overlapind
            case 1
                cross = 1;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_crosshairs
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        crosshair_locations{i}(2*cross-1,seq) = x;
                        crosshair_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 crosshairs
                        dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                        dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <=  90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
            case 2
                cross = 2;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                cross = 1;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 crosshairs
                    dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                    dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
            case 3
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 2nd is now the 3rd cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 3;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                cross = 2;
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                available_points = find(valid_points);
                dist = sqrt((rr-y).^2+(cc-x).^2);
                C = dist <= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                choosen_ind = Cind(randi(length(Cind)));
                [y,x] = ind2sub(size(valid_points),choosen_ind);
                C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                Cind = find(C);
                valid_points(Cind) = 0;
                cross  = 1;%rewrtie to 2 so that when we flip it will work out
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 3:4
                    %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 crosshairs
                    dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                    dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    available_points = find(valid_points);
                    dist = sqrt((rr-y).^2+(cc-x).^2);
                    C = dist <= 15;
                    Cind = find(C);
                    Cind = intersect(Cind,available_points);
                    [Cy,Cx] = ind2sub(size(valid_points),Cind);
                    potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                    for pc = 1:length(Cy)
                        dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                        dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                        angle23 = atan2d(dy23,dx23);
                        potential_angles(pc) = angle23;
                    end
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    good_angles = find(angle_difference <=  90);
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                        crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                        x =  Cx(choosen_ind);
                        y =  Cy(choosen_ind);
                        C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                    end
                end
                xs = crosshair_locations{i}(1:2:end,seq);
                ys = crosshair_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                crosshair_locations{i}(1:2:end,seq) = xs;
                crosshair_locations{i}(2:2:end,seq) = ys;
            case 4
                %going to do this in the reverse order as if this were cross #
                %2, then flip this sequence around so 1st is now the 4th cross.
                %It might be easier to do something else but I already have the
                %code.
                cross = 4;
                x = crosshair_locations{i}(2*cross-1,1);
                y = crosshair_locations{i}(2*cross,1);
                cross = 1;
                crosshair_locations{i}(2*cross-1,seq) = x;
                crosshair_locations{i}(2*cross,seq) = y;
                for cross = 2:number_of_crosshairs
                    if cross == 2
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        choosen_ind = Cind(randi(length(Cind)));
                        [y,x] = ind2sub(size(valid_points),choosen_ind);
                        C = sqrt((rr-y).^2+(cc-x).^2) <= 5;
                        Cind = find(C);
                        valid_points(Cind) = 0;
                        crosshair_locations{i}(2*cross-1,seq) = x;
                        crosshair_locations{i}(2*cross,seq) = y;
                    else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                        %find angle of quickest scan path from previous 2 crosshairs
                        dx12 = crosshair_locations{i}(2*cross-3,seq)-crosshair_locations{i}(2*cross-5,seq);
                        dy12 = crosshair_locations{i}(2*cross-2,seq)-crosshair_locations{i}(2*cross-4,seq);
                        angle12 = atan2d(dy12,dx12);
                        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                        
                        available_points = find(valid_points);
                        dist = sqrt((rr-y).^2+(cc-x).^2);
                        C = dist <= 15;
                        Cind = find(C);
                        Cind = intersect(Cind,available_points);
                        [Cy,Cx] = ind2sub(size(valid_points),Cind);
                        potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                        for pc = 1:length(Cy)
                            dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3,seq);
                            dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2,seq);
                            angle23 = atan2d(dy23,dx23);
                            potential_angles(pc) = angle23;
                        end
                        potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                        angle_difference = potential_angles-angle12;
                        angle_difference = abs(angle_difference);
                        good_angles = find(angle_difference <= 90);
                        if isempty(good_angles);
                            break;
                        else
                            choosen_ind = good_angles(randi(length(good_angles)));
                            crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                            crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                            x =  Cx(choosen_ind);
                            y =  Cy(choosen_ind);
                            C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                            Cind = find(C);
                            valid_points(Cind) = 0;
                        end
                    end
                end
                xs = crosshair_locations{i}(1:2:end,seq);
                ys = crosshair_locations{i}(2:2:end,seq);
                xs = xs(end:-1:1);
                ys = ys(end:-1:1);
                crosshair_locations{i}(1:2:end,seq) = xs;
                crosshair_locations{i}(2:2:end,seq) = ys;
        end
    end
    
    if all(~isnan(crosshair_locations{i}))
        crosshair_locations{i}(1:2:end,:) = crosshair_locations{i}(1:2:end,:)-13;
        crosshair_locations{i}(2:2:end,:) = crosshair_locations{i}(2:2:end,:)-10;
        for seq = 1:size(crosshair_locations{i},2);
            xs =crosshair_locations{i}(1:2:end,seq);
            ys =crosshair_locations{i}(2:2:end,seq);
            
            %double check to make sure distances are good
            d = pdist([xs,ys]);
            if any(d < 5)
                disp('error item locations too close')
                [min(d),seq]
            end
            dx = diff(xs);
            dy = diff(ys);
            if any(sqrt(dx.^2+dy.^2) > 15)
                disp('error locations too far appart')
                crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1;
                if length(crossnum) > 1
                    [crossnum;seq]
                else
                    [crossnum seq]
                end
            end
            if any(abs(xs) > 12) || any(abs(ys) > 9)
                disp('error locations out of bounds')
            end
            
            dy12 = ys(2)-ys(1);
            dx12 = xs(2)-xs(1);
            angle12 = atan2d(dy12,dx12);
            
            %angle from crosshair 2 to 3
            dy23 = ys(3)-ys(2);
            dx23 = xs(3)-xs(2);
            angle23 = atan2d(dy23,dx23);
            
            %angle from crosshair 3 to 4
            dy34 = ys(4)-ys(3);
            dx34 = xs(4)-xs(3);
            angle34 = atan2d(dy34,dx34);
            
            angle12(angle12 < 0) = 360+angle12(angle12 < 0);
            angle23(angle23 < 0) = 360+angle23(angle23 < 0);
            angle34(angle34 < 0) = 360+angle34(angle34 < 0);
            
            dang123 = angle23-angle12;
            dang234 = angle34-angle23;
            dang123(dang123 < 0) =  dang123(dang123 < 0) + 360;
            dang234(dang234 < 0) =  dang123(dang234 < 0) + 360;
            dang123(dang123 > 180) = 360-dang123(dang123 > 180);
            dang234(dang234 > 180) = 360-dang234(dang234 > 180);
            if dang123 > 90 || dang234 > 90
                disp('locations not in a forward direction')
                seq
            end
        end
    end
end

all_nans = [];
for s = 1:length(crosshair_locations);
    if any(any(isnan(crosshair_locations{s}))) 
        all_nans = [all_nans s];
    end
end
crosshair_locations(all_nans) = [];
overlap_index(all_nans) = [];
overlap1 = find(overlap_index == 1);
overlap2 = find(overlap_index == 2);
overlap3 = find(overlap_index == 3);
overlap4 = find(overlap_index == 4);
new_crosshair_locations = [];
new_overlap_index = [];
%ensure that the items that overlaps in a sequences rotates evenly through
%1 to 4
for o = 1:min([length(overlap1),length(overlap2),length(overlap3),length(overlap4)]);
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap1(o))];
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap2(o))];
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap3(o))];
    new_crosshair_locations = [new_crosshair_locations crosshair_locations(overlap4(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap1(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap2(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap3(o))];
    new_overlap_index = [new_overlap_index overlap_index(overlap4(o))];
end
crosshair_locations = new_crosshair_locations; %reassign cuz this is the variable the subsequent code uses
%%
line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------ int1';
line2 =' -4    1      1    0.00    0.00      0   1.00  1.00  0.00                0   0   0 l';
line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50        200 200 200 x';
line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';

clr=['255 255   0 x';
    '  0 255 255 x';
    '255   0 255 x';
    '  0 255   0 x';
    '  0   0   0 x'];

itemnumber = 1:99;%changed from item_letter char('M'+(1:15)-1) SDK 3/5/14
itemtype = [1 2 9 12 12 14]; %rectangle circle elipse triangle hexagon rosshair
int1 = [NaN NaN NaN 3 6 NaN];
rotation = {[0 45 90],0,[0 45 90],[0 30 180 210],0,45};

rand('seed',6996);
for itm = 1:length(itemnumber)
    
    take = randperm(length(itemtype));
    itemtypeorder = itemtype(take);
    if itemnumber(itm) < 10
          set = ['SEQ0' num2str(itemnumber(itm))];
    else
          set = ['SEQ' num2str(itemnumber(itm))];
    end
    fid = fopen([set '.itm'],'w+');
    
    color_order = [clr(randperm(4),:); clr(end,:)];
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    temp_matrix = crosshair_locations{itm}; %first 2 are controlled sequences, 2nd 2 random
    
    %one center clrchng trial for offset correction
    str='  1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  2    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  3    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    
    itmnum = 4;
    for seq = 1:size(temp_matrix,2)
        rot = rotation{take(seq)}(randi(length(rotation{take(seq)})));
        for point = 1:size(temp_matrix,1)/2;
            x_pos(point) = temp_matrix(2*point-1,seq);
            y_pos(point) = temp_matrix(2*point,seq);
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            if itemtypeorder(seq) >= 10
                type_space = '   ';%3 spaces
            else
                type_space = '    ';%3 spaces
            end
            filled_space = '      ';%5 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            if rot < 100
                rotationspace = '  ';
                if rot < 10
                    rotation_add = '.00';
                else
                    rotation_add = '.0';
                end
            else
                rotationspace = '   ';
                rotation_add = '';
            end
            
            if any(itemtypeorder(seq) == [1,9]);
                height_width_space = '0.75  0.38';
            elseif any(itemtypeorder(seq) == [2,12]);
                height_width_space = '0.38  0.38';
            else
                height_width_space = '0.75  0.75';
            end
            
            integer_space = '                        ';
            if any(itemtypeorder(seq) == [1,9,14]);
                inneroutter_space = '              ';
            elseif any(itemtypeorder(seq) == [2,12]);
                inneroutter_space = '  0.38        ';
            else
                inneroutter_space = '  0.75        ';
            end
            if isnan(int1(take(seq)))
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) '\r\n'];
            else
                str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(seq)) filled_space '1' x_space ...
                    num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                    num2str(rot) rotation_add inneroutter_space color_order(seq,:) integer_space num2str(int1(take(seq))) '\r\n'];
            end
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            type_space = '   ';%3 spaces
            filled_space = '      ';%6 spaces
            
            if x_pos(point) <= -10;
                x_space = '  ';%2 spaces
            elseif x_pos(point) < 0
                x_space = '   ';%3 spaces
            elseif x_pos(point) >= 10;
                x_space = '   ';%3 spaces
            else
                x_space = '    ';%4 spaces
            end
            if y_pos(point) <= -10;
                y_space = '  ';%2 spaces
            elseif y_pos(point) < 0
                y_space = '   ';%3 spaces
            elseif y_pos(point) >= 10;
                y_space = '   ';%3 spaces
            else
                y_space = '    ';%4 spaces
            end
            
            str = [itmspace num2str(itmnum) type_space  ' 1' filled_space '1' x_space ...
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   5.00  5.00  0.00'...
                '              ' color_order(end,:) '\r\n'];
            
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
    end
    fclose(fid);
end
%% Create unique random conditions files for each sequence
line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
line2 = '  1     -3    1       1                          2      3'; %color change line

rand('seed',442)
familiarization_block = true;
familiar_trials = 10;
%set to true if you want sequences to be reapeat several times in a row
%during a familrization block
max_conditions = 424;
for itm = 1:length(itemnumber)
    if itemnumber(itm) < 10
        set = ['SEQ0' num2str(itemnumber(itm))];
    else
        set = ['SEQ' num2str(itemnumber(itm))];
    end
    fid = fopen([set '.cnd'],'w+');
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
   
    num_sequences = size(crosshair_locations{itm},2);
    all_conditions = cell(1,num_sequences);
    trials_per_seq = max_conditions/num_sequences;
    
    %for predictable sequnces
    for seq = 1:size(crosshair_locations{itm},2);
        for t = 1:trials_per_seq;
            if seq == 1;
                all_conditions{1,seq} = [all_conditions{1,seq}; 4:11];
            elseif seq == 2;
                all_conditions{1,seq} = [all_conditions{1,seq}; 12:19];
            elseif seq == 3;
                all_conditions{1,seq} = [all_conditions{1,seq}; 20:27];
            elseif seq == 4;
                all_conditions{1,seq} = [all_conditions{1,seq}; 28:35];
            end
        end
    end
    
    %organize conditions as desired
    conditions = [];
    %put into familirization block if desired
    if familiarization_block
        order = randperm(numel(all_conditions));
        for i=1:numel(all_conditions);
            conditions = [conditions;all_conditions{order(i)}(1:familiar_trials,:)];
            all_conditions{order(i)}(1:familiar_trials,:) = [];
        end
    end
    for i=1:numel(all_conditions);
        conditions = [conditions;all_conditions{i}];
    end
    if familiarization_block
        order = 1:size(all_conditions,1)*familiar_trials*num_sequences;
        order = [order randperm(size(conditions,1)-length(order))+length(order)];
        conditions = conditions(order,:);
    else
        conditions = conditions(randperm(size(conditions,1)),:);
    end
    
    %write to file
    for cndline = 1:max_conditions
        if cndline < 9
            cndspace = '  ';
        elseif cndline < 99
            cndspace = ' ';
        else
            cndspace = '';
        end
        btfc = '     -3    2                             '; %background timing, fixid, color palate
        teststr = [];
        for t = 1:8;
            if conditions(cndline,t) < 10;
                testspace = '     ';
            else
                testspace = '    ';
            end
            teststr = [teststr testspace num2str(conditions(cndline,t))];
        end
        fprintf(fid,[cndspace num2str(cndline+1) btfc teststr '\r\n']);
    end
    fclose(fid);
end
