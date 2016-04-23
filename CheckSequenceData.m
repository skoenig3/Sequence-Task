function CheckSequenceData(seq_cortexfile,clrchng_cortexfile,itemnum,numpoints)
% written 6/24/2014 by Seth Konig
% Imports a color change calibration file and cortex file from the sequence
% task into Matlab. Then the function parses the task data (i.e. when cross
% hairs come on and off and detects fixations and saccades using Cluster
% Fix.
%
% Inputs:
%   1) seq_cortexfile: cortex file for the sequence task
%   2) clrchng_cortexfile: cortex file for the clrchng calibration task
%   3) Item number for the sequence task
%   4) numpoints: Number of calibration points (changes which item file to use)
%
% Outputs:
%   1) Matlab file with the save fixations and saccade times as well the
%   sequence task data containing the times at which the crosshairs turned
%   on and off

samprate = 5;

if strcmpi(clrchng_cortexfile(1:2),'PW')
    clrchng_cortexfile = ['R:\Buffalo Lab\Cortex Data\Vivian\' clrchng_cortexfile];
elseif strcmpi(clrchng_cortexfile(1:2),'TT')
    clrchng_cortexfile = ['R:\Buffalo Lab\Cortex Data\Timmy\' clrchng_cortexfile];
end

%%----Color Change Calibration----%%
if numpoints == 25
    ITMFile = 'R:\Buffalo Lab\eblab\Cortex Programs\ClrChng\cch25.itm';
    CNDFile = 'R:\Buffalo Lab\eblab\Cortex Programs\ClrChng\cch25.cnd';
    % this is different becasue the spacing is different and I don't have
    % a new item file on the network for the new spacing
    ind_spacex = [-6,-3,0,3,6]; %whats on the network
    ind_spacey = [-6,-3,0,3,6];%whats on the network
    spacex = [-12,-6,0,6,12];%what actually gets displayed
    spacey = [-8,-4,0,4,8];%what actually gets displayed
    [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = ...
        get_ALLdata(clrchng_cortexfile);
elseif numpoints == 63
    ITMFile = 'C:\Users\seth.koenig\Documents\MATLAB\SSCM\Itm Files\cchgrid.itm';
    CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\SSCM\Cnd Files\cchDrew.cnd';
    ind_spacex = [-16,-12,-8,-4,0,4,8,12,16];
    ind_spacey = [-12,-8,-4,0,4,8,12];
    spacex = [-16,-12,-8,-4,0,4,8,12,16];
    spacey = [-12,-8,-4,0,4,8,12];
    if iscell(clrchng_cortexfile)
        [time_arr1,event_arr1,eog_arr1,epp_arr1,header1,trialcount1]  = ...
            get_ALLdata(['R:\Buffalo Lab\Cortex Data\Vivian\' clrchng_cortexfile{1}]);
        [time_arr2,event_arr2,eog_arr2,epp_arr2,header2,trialcount2]  = ...
            get_ALLdata(['R:\Buffalo Lab\Cortex Data\Vivian\' clrchng_cortexfile{2}]);
        time_arr = [time_arr1 time_arr2];
        event_arr = [event_arr1 event_arr2];
        minsize = min(size(eog_arr1,1),size(eog_arr2,1));
        if size(eog_arr1,1) < size(eog_arr2,1)
            eog_arr1 = [eog_arr1; NaN(size(eog_arr2,1)-size(eog_arr1,1),size(eog_arr1,2))];
        else
            eog_arr2 = [eog_arr2; NaN(size(eog_arr1,1)-size(eog_arr2,1),size(eog_arr2,2))];
        end
        eog_arr = [eog_arr1 eog_arr2];
        trialcount = trialcount1+trialcount2;
    else
        [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(clrchng_cortexfile);
    end
end

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) <= 189
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            if length( perbegind) > 1
                perbegind = perbegind(2);
                perendind = perendind(2);
            end
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                clrchgind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
            end
        end
    end
end

%Don't keep first 2 calibration pionts these are for offset correction at
%start of task
clear cnd
numrpt = size(per,2);
cnd = zeros(1,numrpt);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end

% Create structures x and y of the corresponding average eye data for each trial
% instance (l) of each condition (k)

x = cell(length(spacex),length(spacey));%---For Calibration with Eye tracking data with cp2tform---%
y = cell(length(spacex),length(spacey));
control = NaN(length(cnd),2);
clr = ['rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk'];
figure
hold on
for k = 1:length(cnd)
    C = textscan(itmfil(itmlist(cnd(k)-1000)+5,:),'%d');
    control(k,:) = C{1}(9:10)';
    
    xi = find(C{1}(9) == ind_spacex);
    yi = find(C{1}(10) == ind_spacey);
    eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
    evenind = eyeind(logical(~rem(eyeind,2)));
    oddind =  eyeind(logical(rem(eyeind,2)));
    x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
    y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
end
if iscell(clrchng_cortexfile)
    title(['Calibration transformation for ' clrchng_cortexfile{1}(end-9:end)])
else
    title(['Calibration transformation for ' clrchng_cortexfile(end-9:end)])
end

%Test for errors%
count = zeros(length(spacey),length(spacex));
for xi = 1:length(spacex);
    for yi = 1:length(spacey);
        count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
    end
end
if any(count < 5);
    disp('Calibration trial analysis incomplete or error')
    disp('Check number of calibration pionts or task not finished')
end

clear meanx meany
for k=1:numel(x)
    xss = x{k};
    low = mean(xss)-std(xss);
    high = mean(xss)+std(xss);
    xss(xss < low) = [];
    xss(xss > high) = [];
    meanx(k)=median(xss);
end
for k=1:numel(y)
    yss = y{k};
    low = mean(yss)-std(yss);
    high = mean(yss)+std(yss);
    yss(yss < low) = [];
    yss(yss > high) = [];
    meany(k)=median(y{k});
end

if numpoints == 25
    controlx = [];
    controly = [];
    for i = 1:length(spacex);
        for ii = 1:length(spacey);
            controly = [controly spacey(i)];
            controlx = [controlx spacex(ii)];
        end
    end
else
    controlx = [];
    controly = [];
    for i = 1:length(spacey);
        for ii = 1:length(spacex);
            controlx = [controlx spacex(ii)];
            controly = [controly spacey(i)];
        end
    end
end

if iscell(clrchng_cortexfile)
    tform = cp2tform([controlx' controly'], [meanx' meany'],'polynomial',4);
else
    tform = cp2tform([controlx' controly'], [meanx' meany'],'affine');
end
tform.forward_fcn = tform.inverse_fcn;

figure
hold on
for i = 1:length(controlx);
    plot(controlx(i),controly(i),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
end
if iscell(clrchng_cortexfile)
    title(['Calibration transformation for ' clrchng_cortexfile{1}(end-9:end)])
else
    title(['Calibration transformation for ' clrchng_cortexfile(end-9:end)])
end
xlim([-17.5 17.5])
ylim([-12.5 12.5])

%%---Sequence Task analysis---%%
if strcmpi(seq_cortexfile(1:2),'PW')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Vivian\' seq_cortexfile];
elseif strcmpi(seq_cortexfile(1:2),'TT')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Timmy\' seq_cortexfile];
end

if ischar(itemnum)
    ITMFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Item Files\CKWENZR' itemnum '.itm'];
    CNDFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Item Files\CKWENZR' itemnum '.cnd'];
else
    if itemnum < 10
        ITMFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Item Files\CKWENZ0' num2str(itemnum) '.itm'];
    else
        ITMFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Item Files\CKWENZ' num2str(itemnum) '.itm'];
    end
    CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Item Files\CKWENZ.CND';
end

[time_arr,event_arr,eog_arr,epp_arr, header,trialcount]  = get_ALLdata(cortexfile);

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

if ischar(itemnum)
    itmlist = zeros(size(cndfil,1)-2,4);
    for i = 3:size(cndfil,1)-1;%first cnd is clrchng
        str = textscan(cndfil(i,:),'%d');
        itmlist(i-2,:) = str{1}(4:7);
    end
    crosshair_locations = cell(2,size(itmlist,1)-1); %x in row 1 y in row 2
    for cnd = 1:size(itmlist,1);
        for itm = 1:size(itmlist,2);
            itmline = itmlist(cnd,itm)+6;
            str = textscan(itmfil(itmline,:),'%d');
            crosshair_locations{1,cnd}(1,itm) = str{1}(4); %xpos
            crosshair_locations{1,cnd}(2,itm) = str{1}(5); %ypos
        end
        if double(itemnum) < 69 % RA-RD
            if itmlist(cnd,itm) < 8
                crosshair_locations{2,cnd} = 1; %control sequence 1
            elseif itmlist(cnd,itm) < 16
                crosshair_locations{2,cnd} = 2; %control sequence 2
            elseif itmlist(cnd,itm) < 36
                crosshair_locations{2,cnd} = 3; %random sequence 1
            else
                crosshair_locations{2,cnd} = 4; %random sequence 2
            end
        else
            if itmlist(cnd,itm) < 8
                crosshair_locations{2,cnd} = 1; %control sequence 1
            else
                crosshair_locations{2,cnd} = 2; %random sequence 1
            end
        end
    end
else
    itmlist = zeros(size(cndfil,1)-2,4);
    for i = 2:size(cndfil,1)-1;%last cnd is clrchng
        str = textscan(cndfil(i,:),'%d');
        itmlist(i-1,:) = str{1}(4:7);
    end
    crosshair_locations = cell(1,size(itmlist,1)); %x in row 1 y in row 2
    for cnd = 1:size(itmlist,1);
        for itm = 1:size(itmlist,2);
            itmline = itmlist(cnd,itm)+6;
            str = textscan(itmfil(itmline,:),'%d');
            crosshair_locations{cnd}(1,itm) = str{1}(4); %xpos
            crosshair_locations{cnd}(2,itm) = str{1}(5); %ypos
        end
    end
end

all_trials = [];
numrpt = size(event_arr,2);
new_eog_arr=[];
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 > 1 %1st block is color change always
        if ~isempty(find(event_arr(:,rptlop) == 209)) && ~isempty(find(event_arr(:,rptlop) == 25)) % if a break fixation occured
            perbegind = find(event_arr(:,rptlop) == 100); % eye data starts at 100; might have predictive looking
            perendind = find(event_arr(:,rptlop) == 101); % eye data stops collecting after rewards so can stop here
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop);
            endtimdum = time_arr(perendind(end),rptlop);
            breaktime = time_arr(find(event_arr(:,rptlop) == 209),rptlop)-begtimdum; 
                        
            test0_on = find(event_arr(:,rptlop) == 23);
            test0_off = find(event_arr(:,rptlop) == 24);
            test1_on = find(event_arr(:,rptlop) == 25);
            test1_off = find(event_arr(:,rptlop) == 26);
            test2_on = find(event_arr(:,rptlop) == 27);
            test2_off = find(event_arr(:,rptlop) == 28);
            test3_on = find(event_arr(:,rptlop) == 29);
            test3_off = find(event_arr(:,rptlop) == 30);
            
            %for looking at trials with break fixations 
            test0_on = time_arr(test0_on,rptlop)-begtimdum;
            if isempty(test0_off);
                test0_off = NaN;
            else
                test0_off = time_arr(test0_off,rptlop)-begtimdum;
            end
            if isempty(test1_off);
                test1_off = NaN;
            else
                test1_off = time_arr(test1_off,rptlop)-begtimdum;
            end
            if isempty(test2_off);
                test2_off = NaN;
            else
                test2_off = time_arr(test2_off,rptlop)-begtimdum;
            end
            if isempty(test3_off);
                test3_off = NaN;
            else
                test3_off = time_arr(test3_off,rptlop)-begtimdum;
            end
            if isempty(test1_on);
                test1_on = NaN;
            else
                test1_on = time_arr(test1_on,rptlop)-begtimdum;
            end
            if isempty(test2_on);
                test2_on = NaN; 
            else
                test2_on = time_arr(test2_on,rptlop)-begtimdum;
            end
            if isempty(test3_on);
                test3_on = NaN;
             else
                test3_on = time_arr(test3_on,rptlop)-begtimdum;
            end
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum(1);
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                per(valrptcnt).blk = event_arr(blknumind,rptlop)-500;
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
                per(valrptcnt).test= ...
                    [test0_on test0_off;...
                    test1_on test1_off;...
                    test2_on test2_off;...
                    test3_on test3_off];
                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                all_trials = [all_trials [event_arr(cndnumind,rptlop)-1000; 1]];
                per(valrptcnt).breaktime = breaktime;
            end
        else
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            all_trials = [all_trials [event_arr(cndnumind,rptlop)-1000; 0]];
        end
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
cnd=[];
teststart = [];
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    picstart=1*samprate;
    picend=per(trlop).endsmpind-per(trlop).begsmpind;
    
    eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
    eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
    cnd(trlop)=per(trlop).cnd;
end

%---Recalibrate and automatically scale eye data---%
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x,y] = tformfwd(tform,x,y);
    eyedat{eye} = [x;y];
end
%fixationstats = ClusterFixation_Short(eyedat);
% Note needed to redo cluster Fix for small amounts of data: over lap in clusters in velocity/acceleration state space
% save(['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\' ...
%     seq_cortexfile(1:8) '_' seq_cortexfile(end) 'f-fixation'],'fixationstats','per',...
%     'crosshair_locations','itemnum','all_trials')
save(['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\' ...
    seq_cortexfile(1:8) '_' seq_cortexfile(end) 'f-fixation'],'per',...
    'crosshair_locations','itemnum','all_trials')

end