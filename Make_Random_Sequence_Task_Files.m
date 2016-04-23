% Script Makes Item and condition files for pilotted Random Sequence tasks 
% in which 2 % sequences out of 4 are always presented the same way and the 
% other 2 sequences each have 8 cross hairs that are presented in a random 
% order 4 at a time. Color changes between random and control sequences randomly.

% for ckwenz RA-RD

load('C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Sequence_Task_Sequences.mat');

% find sequences in the 100-109 sequence matrix (since these can't be used
% in the conventional naming of the sequence itm files, and make sure they
% have sequences with a spacing of at least 5 dva between consecuitive
% cross hairs in a given sequence (these are a better distance apart than 3
% dva).
use_these_sequences = [];
for i = 100:length(sequence_matrix);
    valid_positions = ones(19,25);
    good_seq = NaN(1,4);
    for seq = 1:4;
        good_cross = NaN(1,3);
        for cross = 2:4;
            dist_btw_cross = sqrt((sequence_matrix{i}(cross*2-1,seq)-sequence_matrix{i}(cross*2-3,seq)).^2 + ...
                (sequence_matrix{i}(cross*2,seq)-sequence_matrix{i}(cross*2-2,seq)).^2);
            if dist_btw_cross >= 5
                good_cross(cross) = 1;
            else
                good_cross(cross) = 0;
            end
        end
        if all(good_cross);
            good_seq(seq) = 1;
        else
            good_seq(seq) = 0;
        end
    end
    if nansum(good_seq) >= 2
        ind = find(good_seq == 1);
        use_these_sequences = [use_these_sequences; ind(1:2)];
    else
        use_these_sequences = [use_these_sequences; NaN(1,2)];
    end
end
%%
r_sequence_matrix = {};
for i = 1:size(use_these_sequences,1);
    rand('seed',i)
    free_positions = ones(19,25);
    pos = sequence_matrix{i+99}(:,use_these_sequences(i,:));
    r_sequence_matrix{1,i} = sequence_matrix{i+99}(:,use_these_sequences(i,:));
    for seq = 1:2;
        for cross = 1:4
            pos_x = pos(2*cross-1,seq)+13;
            pos_y = pos(2*cross,seq)+10;
            remove_x = pos_x-2:pos_x+2;
            remove_y = pos_y-2:pos_y+2;
            remove_x(remove_x < 1) = [];
            remove_x(remove_x>25) = [];
            remove_y(remove_y < 1) = [];
            remove_y(remove_y>19) = [];
            [x,y] = meshgrid(remove_x,remove_y);
            ind = sub2ind([19,25],y,x);
            free_positions(ind) = 0;
        end
    end
    random_cross = [];%col x,col y, row by cross
    for row = 1:size(free_positions,1);
        if row == 1 || row == size(free_positions,1);
            if rand < 0.5
                continue
            end
        end
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
    else
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
line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------';
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

itemletter = 'ABCDEFGHI';
for itm = 1:9
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
    for sq = 1:size(temp_matrix,2);
        for point = 1:size(temp_matrix,1)/2;
            x_pos(point) = temp_matrix(2*point-1,sq);
            y_pos(point) = temp_matrix(2*point,sq);
        end
        
        for point = 1:length(x_pos);
            if itmnum < 10
                itmspace = '  '; %2 spaces
            else
                itmspace = ' ';%1 space
            end
            type_space = '   ';%3 spaces
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
            
            str = [itmspace num2str(itmnum) type_space  '14' filled_space '1' x_space ...
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   0.75  0.75  45.0'...
                '              ' color_order(sq,:) '\r\n'];
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
        type_space = '   ';%3 spaces
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
        
        if rcross <= 8
            str = [itmspace num2str(itmnum) type_space  '14' filled_space '1' x_space ...
                num2str(x_pos) '.00' y_space num2str(y_pos) '.00      0   0.75  0.75  45.0'...
                '              ' color_order(3,:) '\r\n'];
        else
            str = [itmspace num2str(itmnum) type_space  '14' filled_space '1' x_space ...
                num2str(x_pos) '.00' y_space num2str(y_pos) '.00      0   0.75  0.75  45.0'...
                '              ' color_order(4,:) '\r\n'];
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
%%
% Create Blank item file so can create a save file but it won't run unless
% you load the correct item file for the day

total_items = 51;
fid = fopen('blank.itm','w+');
for line = 1:6
    fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
end
for l = 1:total_items;
    if l < 10
        itmspace = '  '; %2 spaces
    else
        itmspace = ' ';%1 space
    end
    str = [itmspace num2str(l)  '\r\n'];
    fprintf(fid,str);
end
fclose(fid);
%% Create unique random conditions files for each sequence

line1 = 'COND# BCKGND TIMING FIX_ID ---COLOR-PALETTE--- TEST0 TEST1 TEST2 TEST3 TEST4 TEST5 TEST6 TEST7'; %header
line2 = '  1     -3    1       1                          2      3'; %color change line

itemletter = 'ABCDEFGHI';
for itm = 1:9
    set = ['CKWENZR' itemletter(itm)];
    fid = fopen([set '.cnd'],'w+');
    for line = 1:2
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    cndline = 2;
    % "learning" blocks
    seq_order = randperm(4);
    for seq = 1:4;
        for cnd = 1:10
            if cndline < 10
                cndspace = '  ';
            else
                cndspace = ' ';
            end
            btfc = '     -3    2                             '; %background timing, fixid, color palate
            if  seq_order(seq) == 1%control sequences 1
                testline = 4:11; %items and their fixwins
            elseif seq_order(seq) == 2%control sequences 2
                testline = 12:19; %items and their fixwins
            else
                if seq_order(seq) == 3; %random sequence 1;
                    test = 1:8; %20:25
                    r = randperm(8);
                    test = test(r(1:4));
                    testline = [2*test+20-2 test*2+1+20-2];
                elseif seq_order(seq) == 4; %random sequence ;
                    test = 1:8; %starts at 36:51
                    r = randperm(8);
                    test = test(r(1:4));
                    testline = [2*test+36-2 test*2+1+36-2];
                end
            end
            teststr = [];
            for t = 1:8;
                if testline(t) < 10;
                    testspace = '     ';
                else
                    testspace = '    ';
                end
                teststr = [teststr testspace num2str(testline(t))];
            end
            fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
            cndline = cndline+1;
        end
    end
    
    % "testing" blocks
    for cnd = 1:96
        seq_order = randperm(4);
        for seq = 1:4;
            if cndline < 10
                cndspace = '   ';%2 spaces
            elseif cndline < 100
                cndpsace = ' ';% 1 space
            else
                cndspace = '';% 0 space
            end
            btfc = '     -3    2                             '; %background timing, fixid, color palate
            if  seq_order(seq) == 1%control sequences 1
                testline = 4:11; %items and their fixwins
            elseif seq_order(seq) == 2%control sequences 2
                testline = 12:19; %items and their fixwins
            else
                if seq_order(seq) == 3; %random sequence 1;
                    test = 1:8; %20:25
                    r = randperm(8);
                    test = test(r(1:4));
                    testline = [2*test+20-2 test*2+1+20-2];
                elseif seq_order(seq) == 4; %random sequence ;
                    test = 1:8; %starts at 36:51
                    r = randperm(8);
                    test = test(r(1:4));
                    testline = [2*test+36-2 test*2+1+36-2];
                end
            end
            teststr = [];
            for t = 1:8;
                if testline(t) < 10;
                    testspace = '     ';
                else
                    testspace = '    ';
                end
                teststr = [teststr testspace num2str(testline(t))];
            end
            fprintf(fid,[cndspace num2str(cndline) btfc teststr '\r\n']);
            cndline = cndline+1;
        end
    end
    fclose(fid);
end
%%  blank cnd
% Create Blank condition file so can create a save file but it won't run unless
% you load the correct item file for the day

total_cond = 425;
fid = fopen('blank.cnd','w+');
for line = 1
    fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
end
for l = 1:total_cond;
    if l < 10
        cndspace = '   ';%2 spaces
    elseif l < 100
        cndpsace = ' ';% 1 space 
    else
        cndspace = '';% 0 space
    end
    str = [cndspace num2str(l)  '\r\n'];
    fprintf(fid,str);
end
fclose(fid);
