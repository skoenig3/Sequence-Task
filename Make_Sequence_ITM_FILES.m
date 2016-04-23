%% Matlab Code to Automatically Make ITEM files for the Sequence Task
% Generates sequences for ckwenz 5-25 (and 1-4 have to be done manually)
% all sequences are too be random and non overlaping.
% uses a bruit force method since fast to code and fast to computer
% Items 1-4 have to be made manually from item 5 by just creating a single
% sequence in each item file, and each item file is 1 of the 4 sequences in
% item 5.

%%%---Make a bunch of random sequences---%%%
x_pos = -12:12;
y_pos = -9:9;
sequence_length = 4;
num_points = 4;

num_sequences = 10000;
sequence_matrix = cell(1,num_sequences);
for ns = 1:num_sequences;
    temp_matrix = NaN(sequence_length*2,num_points); %x1 y1 x2 y2 ... by row sequence # by col
    for r = 1:size(temp_matrix,1);
        for c = 1:size(temp_matrix,2);
            if rem(r,2) == 1; %odd rows so x pos
                temp_matrix(r,c) = x_pos(randi(length(x_pos),1));
            else %even rows so y pos
                temp_matrix(r,c) = y_pos(randi(length(y_pos),1));
            end
        end
    end
    
    %just write over the positions you want to keep constant
    ox = x_pos(randi(length(x_pos),1));
    oy = y_pos(randi(length(y_pos),1));
    if rem(ns,4) == 0; %overlap position is 1
        temp_matrix(1,1) = ox;
        temp_matrix(2,1) = oy;
        temp_matrix(1,2) = ox;
        temp_matrix(2,2) = oy;
        temp_matrix(1,3) = ox;
        temp_matrix(2,3) = oy;
        temp_matrix(1,4) = ox;
        temp_matrix(2,4) = oy;
    elseif rem(ns,4) == 1; %overlap position is 2
        temp_matrix(3,1) = ox;
        temp_matrix(4,1) = oy;
        temp_matrix(3,2) = ox;
        temp_matrix(4,2) = oy;
        temp_matrix(3,3) = ox;
        temp_matrix(4,3) = oy;
        temp_matrix(3,4) = ox;
        temp_matrix(4,4) = oy;
    elseif rem(ns,4) == 2; %overlap position is 3
        temp_matrix(5,1) = ox;
        temp_matrix(6,1) = oy;
        temp_matrix(5,2) = ox;
        temp_matrix(6,2) = oy;
        temp_matrix(5,3) = ox;
        temp_matrix(6,3) = oy;
        temp_matrix(5,4) = ox;
        temp_matrix(6,4) = oy;
    elseif rem(ns,4) == 3; %overlap position is 4
        temp_matrix(7,1) = ox;
        temp_matrix(8,1) = oy;
        temp_matrix(7,2) = ox;
        temp_matrix(8,2) = oy;
        temp_matrix(7,3) = ox;
        temp_matrix(8,3) = oy;
        temp_matrix(7,4) = ox;
        temp_matrix(8,4) = oy;
    end
    sequence_matrix{ns} = temp_matrix;
end

%%%---Check for good sequences---%%%
%1 done want any points within 3 dva of another point
%2 no overlapping sequences

%check across sequences
temp_matrix = [sequence_matrix{1}(1:2:end)' sequence_matrix{1}(2:2:end)'];
N = size(temp_matrix,1);
[x,y] = meshgrid(1:N);
ind2 = find(ones(N)-eye(N));%find valid pairs not comparisons toself hence -eye the identity matrix
ind2 = [x(ind2) y(ind2)];

bad_sequences = [];
for ns =  1:num_sequences;
    temp_matrix = [sequence_matrix{ns}(1:2:end)' sequence_matrix{ns}(2:2:end)'];
    dist = sqrt((temp_matrix(ind2(:,1),1)-temp_matrix(ind2(:,2),1)).^2+...
        (temp_matrix(ind2(:,1),2)-temp_matrix(ind2(:,2),2)).^2);
    if sum(dist == 0) > 12 %means more than just the builtin point has overlaps
        bad_sequences = [bad_sequences ns];
        continue
    end
    dist(dist == 0) = [];
    if any(dist < 3)
        bad_sequences = [bad_sequences ns];
        continue
    end
end
sequence_matrix(bad_sequences) = [];

%%%---Remove repeat sequences---%%%
count = 1;
while count <= length(sequence_matrix);
    bad_sequences = [];
    for s  = count+1:length(sequence_matrix);
        if all(all(sequence_matrix{count} == sequence_matrix{s}))
            bad_sequences = [bad_sequences s];
        end
    end
    sequence_matrix(bad_sequences) = [];
    count = count+1;
end

%%%---Try to arrange but having dissimilar itm files near each other---%%%
similarity = NaN(length(sequence_matrix),length(sequence_matrix));
count = 1;
while count <= length(sequence_matrix);
    for s  = count+1:length(sequence_matrix);
        similarity(s,count) = sqrt(sum(sum((sequence_matrix{count}-sequence_matrix{s}).^2)));
    end
    count = count+1;
end

%%%---Try to arrange sequences so they're far apart or as dissimilar from day to day---%%%
similarity_thresh = nanmean(nanmean(similarity))+2*nanstd(similarity(1:end));
new_order = 1;
validpoints = 1:length(sequence_matrix);
for s = 1:length(sequence_matrix);
    mx = find(similarity(validpoints,s) == nanmax(similarity(validpoints,s)));
    if ~isempty(mx)
        if length(mx) > 1;
            mx = mx(1);
        end
        mxi = find(similarity(:,s) == similarity(validpoints(mx),s));
        if length(mxi) > 1
            mxi = mxi(1);
        end
        if ~isempty(validpoints)
            new_order = [new_order mxi];
            validpoints(validpoints == mxi) = [];
        end
    end
end
sequence_matrix = sequence_matrix(new_order);

save('Sequence_Task_Sequences.mat','sequence_matrix')
%%
load('C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Sequence_Task_Sequences.mat');

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

num_items = 20;%up to 99
for itm = 5:5+num_items %first 4 are for training
    if itm < 10
        set = ['CKWENZ0' num2str(itm)]; 
    else
        set = ['CKWENZ' num2str(itm)];
    end
    
    fid = fopen([set '.itm'],'w+');
    
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    
    temp_matrix = sequence_matrix{itm-4}; %first 4 are for training
    
    itmnum = 1;
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
                '              ' clr(sq,:) '\r\n'];
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
                num2str(x_pos(point)) '.00' y_space num2str(y_pos(point)) '.00      0   4.00  4.00  0.00'...
                '              ' clr(end,:) '\r\n'];
            fprintf(fid,str);
            itmnum = itmnum+1;
        end
    end
    %one center clrchng trial for offset correction
    str=' 33    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str=' 34    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        150 150 150 x\r\n';
    fprintf(fid,str);
    str=' 35    1      1    0.00    0.00      0   0.30  0.30  0.00  0.30        175 175 130 x\r\n';
    fprintf(fid,str);
    fclose(fid);
end