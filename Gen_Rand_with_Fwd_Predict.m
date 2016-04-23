% Generates sequenes for ckwenz RG - RL
% Written by Seth Konig July 2014
% If trying to generate sequences again RG-RI were generated with an 
% error in which distance between the 2nd and the 3rd, and the 3rd to the 4th 
% items were all referenced to the 2nd item resulting in clustering of the
% 2nd-4th items. These items were still 5-15 dva apart and had a forward
% bias. Sorry. This is what happens when you try to generate new tasks
% files with new parameters really fast. Atleast technically it was still
% producing sequences with the right statistics. 
%---------------------------------------------------------------------%
% to generate G must change line 479 to remove familiratization block %
%---------------------------------------------------------------------%
%load('C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Sequence_Task_Sequences.mat')
Alphabet=char('G'+(1:6)-1);
number_of_sequences = 1;
number_of_crosshairs = 4;
crosshair_locations = {};
buffer = 5; %minimum distance between 2 crosshairs
maxdist = 15; %maximum distance between 2 crosshairs
buffer = buffer - 1;% since items are 1 dva grid
sizex = 25;
sizey = 19;
clear x y
[cc,rr] = meshgrid(1:sizex,1:sizey);
rand('seed',42); %oh yeah! 42!
for i = 1:100; %when this is 100, 42 legitmate sequences actually pop out so awesome
    valid_points = ones(sizey,sizex);
    crosshair_locations{i} = NaN(number_of_crosshairs*2,number_of_sequences);
    for seq = 1:number_of_sequences;
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
                dx12 = crosshair_locations{i}(2*cross-3)-crosshair_locations{i}(2*cross-5);
                dy12 = crosshair_locations{i}(2*cross-2)-crosshair_locations{i}(2*cross-4);
                angle12 = atan2d(dy12,dx12);
                angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                
                available_points = find(valid_points);
                C = sqrt((rr-y).^2+(cc-x).^2)<= 15;
                Cind = find(C);
                Cind = intersect(Cind,available_points);
                [Cy,Cx] = ind2sub(size(valid_points),Cind);
                potential_angles = NaN(1,length(Cy)); %angles between previous crosshair and potential locations
                for pc = 1:length(Cy)
                    dx23 = Cx(pc)-crosshair_locations{i}(2*cross-3);
                    dy23 = Cy(pc)-crosshair_locations{i}(2*cross-2);
                    angle23 = atan2d(dy23,dx23);
                    potential_angles(pc) = angle23;
                end
                potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                angle_difference = potential_angles-angle12;
                angle_difference = abs(angle_difference);
                good_angles = find(angle_difference <=  90);
                if isempty(good_angles);
                    crosshair_locations{i} = NaN(number_of_crosshairs*2,number_of_sequences);
                    break;
                else
                    choosen_ind = good_angles(randi(length(good_angles)));
                    crosshair_locations{i}(2*cross-1,seq) = Cx(choosen_ind);
                    crosshair_locations{i}(2*cross,seq) = Cy(choosen_ind);
                    C = sqrt((rr-Cy(choosen_ind)).^2+(cc-Cx(choosen_ind)).^2) <= 5;
                    
                    %had error in code but didn't notice until vivian had
                    %run through RI.
                    x = Cx(choosen_ind);
                    y = Cy(choosen_ind);
                    Cind = find(C);
                    valid_points(Cind) = 0;
                end
            end
        end
    end
    crosshair_locations{i}(1:2:end,:) = crosshair_locations{i}(1:2:end,:)-13;
    crosshair_locations{i}(2:2:end,:) = crosshair_locations{i}(2:2:end,:)-10;
    xs =crosshair_locations{i}(1:2:end);
    ys =crosshair_locations{i}(2:2:end);
    
    %double check to make sure distances are good
    d = pdist([xs,ys]);
    if any(d < 5)
        disp('error item locations too close')
        min(d)
    end
    dx = diff(xs);
    dy = diff(ys);
    if any(sqrt(dx.^2+dy.^2) > 15)
        disp('error locations too far appart')
        crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1
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
    
    if abs(angle23-angle12) > 90 || abs(angle34-angle23) > 90
        disp('locations not in a forward direction')
    end
end
%%
all_nans = [];
for s = 1:length(crosshair_locations);
    if all(isnan(crosshair_locations{s}))
        all_nans = [all_nans s];
    end
end
crosshair_locations(all_nans) = [];

r_sequence_matrix = [crosshair_locations; cell(1,length(crosshair_locations))];


for i = 1:size(r_sequence_matrix,2);
    rand('seed',42*i); %oh yeah! 42!
    free_positions =ones(sizey,sizex);
    pos = r_sequence_matrix{1,i};
    for seq = 1:size(pos,2)
        for cross = 1:4
            pos_x = pos(2*cross-1,seq)+13;
            pos_y = pos(2*cross,seq)+10;
            C = sqrt((rr-pos_x).^2+(cc-pos_y).^2)<= 5;
            Cind = find(C);
            free_positions(Cind) = 0;
        end
    end
    
    random_cross = [];%col x,col y, row by cross
    for row = 1:size(free_positions,1);
        rline = free_positions(row,:);
        valid_points = find(rline);
        while ~isempty(valid_points);
            random_cross = [random_cross; [valid_points(1),row]];
            remove_x = valid_points(1)-4:valid_points(1)+4;
            remove_y = row-4:row+4;
            remove_x(remove_x < 1) = [];
            remove_x(remove_x>25) = [];
            remove_y(remove_y < 1) = [];
            remove_y(remove_y>19) = [];
            [x,y] = meshgrid(remove_x,remove_y);
            ind = sub2ind([19,25],y,x);
            free_positions(ind) = 0;
            rline = free_positions(row,:);
            valid_points = find(rline);
        end
    end
    random_cross(:,1) = random_cross(:,1) - 13;
    random_cross(:,2) = random_cross(:,2) - 10;
    if size(random_cross,1) >= 16
        take = randperm(size(random_cross,1));
        take = take(1:16);
        random_cross = random_cross(take,:);
    elseif size(random_cross,1) < 8
        disp('Note enough random crosses')
    end
    r_sequence_matrix{2,i} = random_cross;
end

for seq = 1:size(r_sequence_matrix,2)
    if any(abs(r_sequence_matrix{2,seq}(:,1)) > 12 | abs(r_sequence_matrix{2,seq}(:,2)) > 9)
        disp([num2str(seq) ' does not have a valid point. Double check manually'])
    end
end
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

itemletter = char('G'+(1:6)-1);
itemtype = [1 2 9 12 12 14]; %rectangle circle elipse triangle hexagon rosshair
int1 = [NaN NaN NaN 3 6 NaN];
rotation = {[0 45 90],0,[0 45 90],[0 30 180 210],0,45};

rand('seed',10101);
for itm = 1:6%1:9
    
    take = randperm(length(itemtype));
    take = take(1:2);
    itemtypeorder = itemtype(take);
    rot1 = rotation{take(1)}(randi(length(rotation{take(1)})));
    rot2 = rotation{take(2)}(randi(length(rotation{take(2)})));
    set = ['CKWENZR' itemletter(itm)];
    fid = fopen([set '.itm'],'w+');
    
    color_order = [clr(randperm(4),:); clr(end,:)];
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    temp_matrix = r_sequence_matrix{1,itm}; %first 2 are controlled sequences, 2nd 2 random
    
    %one center clrchng trial for offset correction
    str='  1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  2    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str='  3    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    
    itmnum = 4;
    for point = 1:size(temp_matrix,1)/2;
        x_pos(point) = temp_matrix(2*point-1,1);
        y_pos(point) = temp_matrix(2*point,1);
    end
    
    for point = 1:length(x_pos);
        if itmnum < 10
            itmspace = '  '; %2 spaces
        else
            itmspace = ' ';%1 space
        end
        if itemtypeorder(1) >= 10
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
        
        if rot1 < 100
            rotationspace = '  ';
            if rot1 < 10
                rotation_add = '.00';
            else
                rotation_add = '.0';
            end
        else
            rotationspace = '   ';
            rotation_add = '';
        end
        
        if any(itemtypeorder(1) == [1,9]);
            height_width_space = '0.75  0.38';
        elseif any(itemtypeorder(1) == [2,12]);
            height_width_space = '0.38  0.38';
        else
            height_width_space = '0.75  0.75';
        end
        
        integer_space = '                        ';
        if any(itemtypeorder(1) == [1,9,14]);
            inneroutter_space = '              ';
        elseif any(itemtypeorder(1) == [2,12]);
            inneroutter_space = '  0.38        ';
        else
            inneroutter_space = '  0.75        ';
        end
        if isnan(int1(take(1)))
            str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(1)) filled_space '1' x_space ...
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                num2str(rot1) rotation_add inneroutter_space color_order(1,:) '\r\n'];
        else
            str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(1)) filled_space '1' x_space ...
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   ' height_width_space rotationspace...
                num2str(rot1) rotation_add inneroutter_space color_order(1,:) integer_space num2str(int1(take(1))) '\r\n'];
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
    
    %for the random cross hairs
    
    temp_matrix = r_sequence_matrix{2,itm}; %first 2 are controlled sequences, 2nd 2 random
    for rcross = 1:size(temp_matrix,1);
        x_pos = temp_matrix(rcross,1);
        y_pos = temp_matrix(rcross,2);
        
        if itmnum < 10
            itmspace = '  '; %2 spaces
        else
            itmspace = ' ';%1 space
        end
        if itemtypeorder(2) >= 10
            type_space = '   ';%3 spaces
        else
            type_space = '    ';%3 spaces
        end
        
        filled_space = '      ';%5 spaces
        
        if x_pos <= -10;
            x_space = '  ';%2 spaces
        elseif x_pos < 0
            x_space = '   ';%3 spaces
        elseif x_pos >= 10;
            x_space = '   ';%3 spaces
        else
            x_space = '    ';%4 spaces
        end
        if y_pos <= -10;
            y_space = '  ';%2 spaces
        elseif y_pos < 0
            y_space = '   ';%3 spaces
        elseif y_pos >= 10;
            y_space = '   ';%3 spaces
        else
            y_space = '    ';%4 spaces
        end
        
        if any(itemtypeorder(2) == [1,9]);
            height_width_space = '0.75  0.38';
        elseif any(itemtypeorder(2) == [2,12]);
            height_width_space = '0.38  0.38';
        else
            height_width_space = '0.75  0.75';
        end
        
        if rot2 < 100
            rotationspace = '  ';
            if rot2 < 10
                rotation_add = '.00';
            else
                rotation_add = '.0';
            end
        else
            rotationspace = '   ';
            rotation_add = '';
        end
        
        integer_space = '                        ';
        if any(itemtypeorder(2) == [1,9,14]);
            inneroutter_space = '              ';
        elseif any(itemtypeorder(2) == [2,12]);
            inneroutter_space = '  0.38        ';
        else
            inneroutter_space = '  0.75        ';
        end
        if isnan(int1(take(2)))
            str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(2)) filled_space '1' x_space ...
                num2str(x_pos) '.00' y_space num2str(y_pos) '.00      0   '  height_width_space rotationspace...
                num2str(rot2) rotation_add inneroutter_space color_order(2,:) '\r\n'];
        else
            str = [itmspace num2str(itmnum) type_space  num2str(itemtypeorder(2)) filled_space '1' x_space ...
                num2str(x_pos) '.00' y_space num2str(y_pos) '.00      0   ' height_width_space rotationspace...
                num2str(rot2) rotation_add inneroutter_space color_order(2,:) integer_space num2str(int1(take(2))) '\r\n'];
        end
        
        fprintf(fid,str);
        itmnum = itmnum+1;
        
        if itmnum < 10
            itmspace = '  '; %2 spaces
        else
            itmspace = ' ';%1 space
        end
        type_space = '   ';%3 spaces
        filled_space = '      ';%6 spaces
        
        if x_pos <= -10;
            x_space = '  ';%2 spaces
        elseif x_pos < 0
            x_space = '   ';%3 spaces
        elseif x_pos >= 10;
            x_space = '   ';%3 spaces
        else
            x_space = '    ';%4 spaces
        end
        if y_pos <= -10;
            y_space = '  ';%2 spaces
        elseif y_pos < 0
            y_space = '   ';%3 spaces
        elseif y_pos >= 10;
            y_space = '   ';%3 spaces
        else
            y_space = '    ';%4 spaces
        end
        
        str = [itmspace num2str(itmnum) type_space  ' 1' filled_space '1' x_space ...
            num2str(x_pos) '.00' y_space num2str(y_pos) '.00      0   5.00  5.00  0.00'...
            '              ' clr(end,:) '\r\n'];
        fprintf(fid,str);
        itmnum = itmnum+1;
    end
    fclose(fid);
end
%% Create unique random conditions files for each sequence
line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
line2 = '  1     -3    1       1                          2      3'; %color change line

rand('seed',142)
itemletter = char('G'+(1:6)-1);
familiarization_block = true; %false for G
familiar_trials = 10;
%set to true if you want sequences to be reapeat several times in a row
%during a familrization block
max_conditions = 424;
for itm = 1:6
    set = ['CKWENZR' itemletter(itm)];
    fid = fopen([set '.cnd'],'w+');
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    all_conditions = cell(2,size(r_sequence_matrix{1,itm},2));
    num_sequences = size(r_sequence_matrix{1,itm},2)+size(r_sequence_matrix{2,itm},2)/2;
    trials_per_seq = max_conditions/num_sequences;
    
    %for predictable sequnces
    for seq = 1:size(r_sequence_matrix{1,itm},2);
        for t = 1:trials_per_seq;
            if seq == 1;
                all_conditions{1,seq} = [all_conditions{1,seq}; 4:11];
            elseif seq == 2;
                all_conditions{1,seq} = [all_conditions{1,seq}; 12:16];
            elseif seq == 3;
                all_conditions{1,seq} = [all_conditions{1,seq}; 17:24];
            elseif seq == 4;
                all_conditions{1,seq} = [all_conditions{1,seq}; 25:32];
            end
        end
    end
    
    %for random sequences
    if ~isempty(r_sequence_matrix{2,itm}) %if there actually are
        count = 1; %don't know how many tries it's going to take to find valid sequences
        while count <= trials_per_seq;
            locs = r_sequence_matrix{2,itm};
            all_locs = locs;
            crossnum = 1:length(locs);
            test = NaN(1,4);
            for cross = 1:4;
                if cross == 1
                    r = randi(length(locs));
                    test(cross) = r;
                    crossnum(r)  = [];
                    locs(r,:) = [];
                elseif cross == 2
                    dists = sqrt((locs(:,1)-all_locs(test(cross-1),1)).^2 + (locs(:,2)-all_locs(test(cross-1),2)).^2);
                    valid_ind = find(dists <= 15);
                    choosen_ind = valid_ind(randi(length(valid_ind)));
                    locs(choosen_ind,:) = [];
                    test(cross) = crossnum(choosen_ind);
                    crossnum(choosen_ind) = [];
                else %for cross 3 and 4 ensure forward/saccadic momentum doesn't interfere with reaction times
                    %find angle of quickest scan path from previous 2 crosshairs
                    dists = sqrt((locs(:,1)-all_locs(test(cross-1),1)).^2 + (locs(:,2)-all_locs(test(cross-1),2)).^2);
                    valid_ind = find(dists <= 15);
                   
                    dx12 = all_locs(test(cross-1),1)-all_locs(test(cross-2),1);
                    dy12 = all_locs(test(cross-1),2)-all_locs(test(cross-2),2);
                    angle12 = atan2d(dy12,dx12);
                    angle12(angle12 < 0) = 360+angle12(angle12 < 0);
                    
                    dx23 = locs(valid_ind,1)-all_locs(test(cross-1),1);
                    dy23 = locs(valid_ind,2)-all_locs(test(cross-1),2);
                    potential_angles  = atan2d(dy23,dx23);
                    potential_angles(potential_angles < 0) = 360+potential_angles(potential_angles < 0);
                    angle_difference = potential_angles-angle12;
                    angle_difference = abs(angle_difference);
                    angle_difference(angle_difference > 180) = 360-angle_difference(angle_difference > 180);
                    good_angles = find(abs(angle_difference) <=  90);
                    
                    if isempty(good_angles);
                        break;
                    else
                        choosen_ind = good_angles(randi(length(good_angles)));
                        test(cross) = crossnum(valid_ind(choosen_ind));
                        if cross == 3
                            crossnum(valid_ind(choosen_ind)) = [];
                            locs(valid_ind(choosen_ind),:) = [];
                        end
                    end
                end
                if cross == 4
                    xs = all_locs(test,1);
                    ys = all_locs(test,2);
                    d = pdist([xs,ys]);
                    if any(d < 5)
                        disp('error item locations too close')
                        min(d)
                    end
                    dx = diff(xs);
                    dy = diff(ys);
                    if any(sqrt(dx.^2+dy.^2) > 15)
                        disp('error locations too far appart')
                        crossnum = find(sqrt(dx.^2+dy.^2) > 15)+1
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
                    
                    d1 = abs(angle23-angle12);
                    d1(d1 > 180) = 360-d1(d1 > 180);
                    d2 = abs(angle34-angle23);
                    d2(d2 > 180) = 360-d2(d2 > 180);
                    if d1 > 90 || d1 > 90
                        disp('locations not in a forward direction')
                        d1
                        d2
                    end
                    all_conditions{2,1} = [all_conditions{2,1}; 2*test+12-2 test*2+1+12-2];
                    count = count+1;
                end
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
        order = 1:size(all_conditions,1)*familiar_trials;
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
