% Random Sequence Task Analysis
% Written in June 2014 by Seth König adapted for the random July 1, 2014
% Essentially the same as the "Sequence_Task_Analysis.m" except adatped for
% analysis of the random sequence task pilot study. 
%
% For sequences RA-RD
%
% [1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix
% [2] Combine Behavior Across several different item sets

%%
%%---[1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix---%%
sequence_cortex_files = {'PW140701.4','PW140702.2','PW140702.4','PW140703.2','PW140703.4'};
calibration_cortex_files = {{'PW140701.1','PW140701.3'},'PW140702.1','PW140702.1',...
    'PW140703.1','PW140703.1','PW140707.1'};
item_letter = ['A','D','D','B','B','C'];
number_calibration_points = [63 25 25 25 25 25];

% sequence_cortex_files = {'TT140908.2','TT140910.1','TT140911.2','TT140912.2'};
% calibration_cortex_files =  {'TT140908.1','TT140910.2','TT140911.1','TT140912.1'};
% item_letter = ['ABCD'];
% number_calibration_points = [26 26 26 26];

% sequence_cortex_files = {'RR141024.2','RR141027.2','RR141028.2','RR141029.2'};
% calibration_cortex_files = {'RR141024.1','RR141027.1','RR141028.1','RR141029.1'};
% item_letter = ['ABCD'];
% number_calibration_points = [25*ones(1,4)];

% for file = 1:length(sequence_cortex_files);
%     getSequenceData(sequence_cortex_files{file},calibration_cortex_files{file},...
%         item_letter(file),number_calibration_points(file));
%     %close all
% end
%%
%%---[2] Combine Behavior Across several different item sets---%%
eyedat_dir = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\'];
eye_tracking_error = 5; %dva of fixation window
clr = ['bycy'];

fixation_files = sequence_cortex_files;

condition = cell(1,length(fixation_files));
reaction_time = cell(1,length(fixation_files)); %how long after a stimulus disappears do they saccade
time_to_fixation = cell(1,length(fixation_files)); %how fast did they get to it the first time around
fixation_accuracy = cell(1,length(fixation_files)); %how far off
fixation_duration = cell(1,length(fixation_files)); %fixation duration
fixation_on_test = cell(1,length(fixation_files)); %how long the stimulus was on while fixating
extrafixations = cell(1,length(fixation_files)); %how many noncross hair fixations do they make
for file = 1:length(fixation_files);
    load([eyedat_dir fixation_files{file}(1:8) '_' fixation_files{file}(end) '-fixation.mat'])
    condition{file} = NaN(length(fixationstats),1); %which cross hair was presented
    reaction_time{file} = NaN(length(fixationstats),size(per(1).test,1)-1); %how long after a stimulus disappears do they saccade
    time_to_fixation{file} = NaN(length(fixationstats),size(per(1).test,1)); %how fast did they get to it the first time around
    fixation_accuracy{file} = NaN(length(fixationstats),size(per(1).test,1)); %how far off
    fixation_duration{file} = NaN(length(fixationstats),size(per(1).test,1)); %fixation duration
    fixation_on_test{file} = NaN(length(fixationstats),size(per(1).test,1)); %how long the stimulus was on while fixating
    extrafixations{file} = zeros(length(fixationstats),size(per(1).test,1)); %how many noncross hair fixations do they make
    
    for trial = 1:length(fixationstats)
        if 1%length(find(per(trial).cnd == all_trials(1,:))) == 1%otherwise saw sequence multiple times
            fixations = fixationstats{trial}.fixations;
            fixationtimes = fixationstats{trial}.fixationtimes;
            saccadetimes = fixationstats{trial}.saccadetimes;
            xy = fixationstats{trial}.XY;
            
            % Find saccades that were small in amplitude and determine if these
            % smaller saccades were corrective and I should essentially count
            % the pre and post saccade fixations as both on the cross hair.
            saccade_amplitudes = NaN(1,size(saccadetimes,2));
            for sac = 1:size(saccadetimes,2);
                sacx = xy(1,saccadetimes(2,sac))-xy(1,saccadetimes(1,sac));
                sacy = xy(2,saccadetimes(2,sac))-xy(2,saccadetimes(1,sac));
                saccade_amplitudes(sac) = sqrt(sacx.^2+sacy.^2);
            end
            if any(saccade_amplitudes(2:end) <= eye_tracking_error/2) %going to generalize  these as corrective saccades. Ignore 1st saccade
                [fixationtimes,fixations] = remove_corrective_saccades(...
                    xy,saccadetimes,saccade_amplitudes,eye_tracking_error);
            end
            
            crosses = double(crosshair_locations{1,per(trial).cnd-1}(:,1:size(per(trial).test,1)));
            
            %         figure
            %         hold on
            %         plot(xy(1,:),xy(2,:),'g');
            % %         for f = 1:length(fixationtimes);
            % %             plot(fixations(1,f),fixations(2,f),'k*','markersize',5);
            % %         end
            %         for c = 1:size(crosses,2);
            %             plot(crosses(1,c),crosses(2,c),'kx','markersize',12)
            %         end
            %         xlim([-17 17])
            %         ylim([-14 14])
            
            valid_fixations = ones(1,size(fixations,2));
            for c = 1:size(crosses,2)
                if c == 1
                    valid_window = [per(trial).test(c,1) per(trial).test(c,2)];
                else
                    valid_window = [per(trial).test(c-1,2) per(trial).test(c,2)];
                end
                valid_ind = valid_window(1):valid_window(2);
                valid_ind(valid_ind > length(xy)) = [];
                
                potential_fix = [];
                for f = 1:size(fixationtimes,2);
                    if valid_fixations(f)
                        fixtimes = fixationtimes(1,f):fixationtimes(2,f);
                        C = intersect(fixtimes,valid_ind);
                        if length(C) >= 50;
                            potential_fix = [potential_fix f];
                        end
                    end
                end
                if ~isempty(potential_fix)
                    if length(potential_fix)  == 1;
                        fixstart = fixationtimes(1,potential_fix);
                        valid_fixations(1:potential_fix) = 0;
                        dist = sqrt((fixations(1,potential_fix)-crosses(1,c)).^2+...
                            (fixations(2,potential_fix)-crosses(2,c)).^2);
                    elseif length(potential_fix) > 1;
                        dist = sqrt((fixations(1,potential_fix)-crosses(1,c)).^2+...
                            (fixations(2,potential_fix)-crosses(2,c)).^2);
                        [~,thefix] = min(dist);
                        fixstart = fixationtimes(1,potential_fix(thefix));
                        valid_fixations(1:potential_fix(thefix)) = 0;
                        potential_fix = potential_fix(thefix);
                        dist = dist(thefix);
                        if thefix ~= 1
                            extrafixations{file}(trial,c) = thefix-1;
                        end
                    end
                    
                    if crosshair_locations{2,per(trial).cnd-1} <= 2
                        rand_condition{file}(trial) = 1; %control
                    else
                        rand_condition{file}(trial) = 2;% random
                    end
                    condition{file}(trial)= crosshair_locations{2,per(trial).cnd-1};
                    
                    time_to_fixation{file}(trial,c) = fixstart-per(trial).test(c,1);
                    
                    fixation_duration{file}(trial,c) = fixationtimes(2,potential_fix)-fixationtimes(1,potential_fix)+1;
                    if c < size(crosses,2)
                        reaction_time{file}(trial,c) = fixationtimes(2,potential_fix)-per(trial).test(c,2);
                    end
                    
                    timeON = per(trial).test(c,1):per(trial).test(c,2);
                    fixON = fixationtimes(1,potential_fix):fixationtimes(2,potential_fix);
                    fixation_on_test{file}(trial,c) = length(intersect(timeON,fixON));
                    
                    fixation_accuracy{file}(trial,c) = dist;
                    %                 plot(xy(1,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),...
                    %                     xy(2,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),'c');
                    %                 plot(xy(1,fixstart),xy(2,fixstart),'m*','markersize',10)
                end
            end
            %         pause(0.5)
            %         close
        end
    end
end

all_reaction_time = cell(2,250); %how long after a stimulus disappears do they saccade
all_time_to_fixation = cell(2,250); %how fast did they get to it the first time around
all_fixation_accuracy = cell(2,250); %how far off
all_fixation_duration = cell(2,250); %fixation duration
all_fixation_on_test = cell(2,250); %how long the stimulus was on while fixating
all_extrafixations = cell(2,250); %how many noncross hair fixations do they make
for file = 1:length(condition);
    for cnd = 1:4
        trials = find(condition{file} == cnd);
        for repitition = 1:length(trials)
            if rand_condition{file}(trials(repitition)) == 1;
                rand_row = 1;
            else
                rand_row = 2;
            end
            all_reaction_time{rand_row,repitition} = [all_reaction_time{rand_row,repitition}; reaction_time{file}(trials(repitition),:)];
            all_time_to_fixation{rand_row,repitition} = [all_time_to_fixation{rand_row,repitition}; time_to_fixation{file}(trials(repitition),:)];
            all_fixation_accuracy{rand_row,repitition} = [all_fixation_accuracy{rand_row,repitition}; fixation_accuracy{file}(trials(repitition),:)];
            all_fixation_duration{rand_row,repitition} = [all_fixation_duration{rand_row,repitition}; fixation_duration{file}(trials(repitition),:)];
            all_fixation_on_test{rand_row,repitition} = [all_fixation_on_test{rand_row,repitition}; fixation_on_test{file}(trials(repitition),:)];
            all_extrafixations{rand_row,repitition} = [all_extrafixations{rand_row,repitition}; extrafixations{file}(trials(repitition),:)];
        end
    end
end

too_few = cellfun(@numel,all_reaction_time);
too_few = find(too_few(1,:) <= 0);
all_reaction_time(:,too_few) =[];
all_time_to_fixation(:,too_few)=[];
all_fixation_accuracy(:,too_few)=[];
all_fixation_duration(:,too_few)=[];
all_fixation_on_test(:,too_few)=[];
all_extrafixations(:,too_few)=[];

mean_reaction_time = cell(1,2);
mean_time_to_fixation =cell(1,2);
mean_fixation_accuracy =cell(1,2);
mean_fixation_duration =cell(1,2);
mean_fixation_on_test =cell(1,2);
mean_extrafixations =cell(1,2);
std_reaction_time =cell(1,2);
std_time_to_fixation =cell(1,2);
std_fixation_accuracy =cell(1,2);
std_fixation_duration =cell(1,2);
std_fixation_on_test =cell(1,2);
std_extrafixations =cell(1,2);
num_reaction_time =cell(1,2);
num_time_to_fixation =cell(1,2);
num_fixation_accuracy =cell(1,2);
num_fixation_duration =cell(1,2);
num_fixation_on_test =cell(1,2);
num_extrafixations =cell(1,2);
for lt = 1:length(all_reaction_time)-1;
    for cnd = 1:2 %row 1 is control, row 2 is random
        for cross = 1:4;
            mean_time_to_fixation{cnd}(cross,lt) = nanmean(all_time_to_fixation{cnd,lt}(:,cross));
            std_time_to_fixation{cnd}(cross,lt) = nanstd(all_time_to_fixation{cnd,lt}(:,cross));
            num_time_to_fixation{cnd}(cross,lt) = sum(~isnan(all_time_to_fixation{cnd,lt}(:,cross)));
            mean_fixation_accuracy{cnd}(cross,lt) = nanmean(all_fixation_accuracy{cnd,lt}(:,cross));
            std_fixation_accuracy{cnd}(cross,lt) = nanstd(all_fixation_accuracy{cnd,lt}(:,cross));
            num_fixation_accuracy{cnd}(cross,lt) = sum(~isnan(all_fixation_accuracy{cnd,lt}(:,cross)));
            mean_fixation_on_test{cnd}(cross,lt) = nanmean(all_fixation_on_test{cnd,lt}(:,cross));
            std_fixation_on_test{cnd}(cross,lt) = nanstd(all_fixation_on_test{cnd,lt}(:,cross));
            num_fixation_on_test{cnd}(cross,lt) = sum(~isnan(all_fixation_on_test{cnd,lt}(:,cross)));
            mean_extrafixations{cnd}(cross,lt) = nanmean(all_extrafixations{cnd,lt}(:,cross));
            std_extrafixations{cnd}(cross,lt) = nanstd(all_extrafixations{cnd,lt}(:,cross));
            num_extrafixations{cnd}(cross,lt) = sum(~isnan(all_extrafixations{cnd,lt}(:,cross)));
            mean_fixation_duration{cnd}(cross,lt) = nanmean(all_fixation_duration{cnd,lt}(:,cross));
            std_fixation_duration{cnd}(cross,lt) = nanstd(all_fixation_duration{cnd,lt}(:,cross));
            num_fixation_duration{cnd}(cross,lt) = sum(~isnan(all_fixation_duration{cnd,lt}(:,cross)));
            if cross <= 3
                mean_reaction_time{cnd}(cross,lt) = nanmean(all_reaction_time{cnd,lt}(:,cross));
                std_reaction_time{cnd}(cross,lt) = nanstd(all_reaction_time{cnd,lt}(:,cross));
                num_reaction_time{cnd}(cross,lt) = sum(~isnan(all_reaction_time{cnd,lt}(:,cross)));
            end
        end
    end
end

figure
for cross = 1:4
    subplot(2,2,cross)
    hold on
    errorbar(mean_time_to_fixation{1}(cross,:),std_time_to_fixation{1}(cross,:)./sqrt(num_time_to_fixation{1}(1,:)))
    errorbar(mean_time_to_fixation{2}(cross,:),std_time_to_fixation{2}(cross,:)./sqrt(num_time_to_fixation{2}(1,:)),'r')
    hold off
    xlabel('Trial #')
    ylabel('Time to Fixation (ms)')
    title(['Cross ' num2str(cross)])
    legend('Predictable','Random')
end
subtitle('Time to Fixation')

frt = cell(1,2);
for cnd = 1:2
    for trial = 11:length(all_time_to_fixation);
        frt{cnd} = [frt{cnd};all_time_to_fixation{cnd,trial}];
    end
end
means = NaN(2,4);
stds = NaN(2,4);
nmp = NaN(2,4);
for cnd = 1:2
    for cross = 1:4
        means(cnd,cross) = nanmean(frt{cnd}(:,cross));
        stds(cnd,cross) = nanstd(frt{cnd}(:,cross));
        nmp(cnd,cross) = sqrt(sum(~isnan(frt{cnd}(:,cross))));
    end
end
figure
errorb(means',(stds./nmp)')
ylabel('Time to Fixation (ms)')
xlabel('Cross #')
title('Mean Time to Fixation')
legend('Predictable','Random')

nmp = frt{1}(:,2:end);
nmp = nmp(1:end);
nmp = sqrt(sum(~isnan(nmp)));
figure
hold on
bar(1,nanmean(sum(frt{1}(:,2:end),2)))
bar(2,nanmean(sum(frt{2}(:,2:end),2)),'r')
errorb(1,nanmean(sum(frt{1}(:,2:end),2)),nanstd(sum(frt{1}(:,2:end),2))/nmp)
errorb(2,nanmean(sum(frt{2}(:,2:end),2)),nanstd(sum(frt{2}(:,2:end),2))/nmp)
ylabel('Cummulative Time to Fixation (ms)')
set(gca,'Xtick',[1:2])
set(gca,'XtickLabel',{'Predictable','Random'})
ylim([500 600])

figure
for cross = 1:4
    subplot(2,2,cross)
    hold on
    errorbar(mean_fixation_accuracy{1}(cross,:),std_fixation_accuracy{1}(cross,:)./sqrt(num_fixation_accuracy{1}(1,:)))
    errorbar(mean_fixation_accuracy{2}(cross,:),std_fixation_accuracy{2}(cross,:)./sqrt(num_fixation_accuracy{2}(1,:)),'r')
    hold off
    xlabel('Trial #')
    ylabel('Fixation Accuracy (dva)')
    title(['Cross ' num2str(cross)])
    legend('Predictable','Random')
end
subtitle('Fixation Accuracy')

figure
for cross = 1:4
    subplot(2,2,cross)
    hold on
    errorbar(mean_extrafixations{1}(cross,:),std_extrafixations{1}(cross,:)./sqrt(num_extrafixations{1}(1,:)))
    errorbar(mean_extrafixations{2}(cross,:),std_extrafixations{2}(cross,:)./sqrt(num_extrafixations{2}(1,:)),'r')
    hold off
    xlabel('Trial #')
    ylabel('# of Extra Fixations')
    title(['Cross ' num2str(cross)])
    legend('Predictable','Random')
    ylim([0 2])
end
subtitle('# of Extra Fixations')
%%
%plot individual day data
plots = {1,2,4,5,[3 6]};
for file = 1:length(time_to_fixation);
    figure
    trials = find(rand_condition{file} == 1);
    rtrials = find(rand_condition{file} == 2);
    for cnd = 1:4
        cnd_data = time_to_fixation{file}(trials,cnd);
        rcnd_data = time_to_fixation{file}(rtrials,cnd);
        subplot(2,3,plots{cnd})
        hold on
        plot(cnd_data);
        plot(rcnd_data,'r');
        ylim([-250 750])
        xlabel('Trial #')
        ylabel('Time to Fixation')
        legend('Predictable','Random');
        title(['Cross #' num2str(cnd)])
    end
    means = [nanmean(time_to_fixation{file}(trials,:));nanmean(time_to_fixation{file}(rtrials,:))];
    stds = [nanstd(time_to_fixation{file}(trials,:));nanstd(time_to_fixation{file}(rtrials,:))];
    subplot(2,3,plots{end})
    hold on
    bar(means')
    errorb(means',stds')
    hold off
    ylim([0 350])
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Cross 1','Cross 2','Cross 3','Cross 4'})
    ylabel('Time to Fixation (ms)')
    subtitle(fixation_files{file}(1:8))
end
%% is saccadic momentum involved in reaction times for 3rd and 4th crosshairs?

reaction_time = cell(1,2);
ang = cell(1,2);
for file = 1:length(fixation_files);
    load([eyedat_dir fixation_files{file}(1:8) '_' fixation_files{file}(end) '-fixation.mat'])
    for trial = 1:size(time_to_fixation{file},1);
        if crosshair_locations{2,trial} <= 2
            rand_row = 1; %predictable sequence
        else
            rand_row = 2; %random sequence
        end
        crosses = double(crosshair_locations{1,trial});
        
        %angle from crosshair 1 to 2
        dy12 = crosses(2,2)-crosses(2,1);
        dx12 = crosses(1,2)-crosses(1,1);
        angle12 = atan2d(dy12,dx12);
        
        %angle from crosshair 2 to 3
        dy23 = crosses(2,3)-crosses(2,2);
        dx23 = crosses(1,3)-crosses(1,2);
        angle23 = atan2d(dy23,dx23);
              
        %angle from crosshair 3 to 4
        dy34 = crosses(2,4)-crosses(2,3);
        dx34 = crosses(1,4)-crosses(1,3);
        angle34 = atan2d(dy34,dx34);

        angle12(angle12 < 0) = 360+angle12(angle12 < 0);
        angle23(angle23 < 0) = 360+angle23(angle23 < 0);
        angle34(angle34 < 0) = 360+angle34(angle34 < 0);


        ang{rand_row} = [ang{rand_row} angle23-angle12];
        reaction_time{rand_row} = [reaction_time{rand_row} time_to_fixation{file}(trial,3)];
        ang{rand_row} = [ang{rand_row} angle34-angle23];
        reaction_time{rand_row} = [reaction_time{rand_row} time_to_fixation{file}(trial,4)];
    end
end

ang{1} = abs(ang{1});
ang{1}(ang{1} > 180) = ang{1}(ang{1} > 180)-180;
ang{2} = abs(ang{2});
ang{2}(ang{2} > 180) = ang{2}(ang{2} > 180)-180;

labels = {'10','20','30','40','50','60','70','80','90','100','110','120','130','140','150','160','170','180'}; 

figure
hold on
plot(ang{1},reaction_time{1},'.');
plot(ang{2},reaction_time{2},'r.');
xlabel('Change in Saccade Angle')
ylabel('Time to Fixation (ms)')
legend('Predictable','Random')
set(gca,'Xtick',0:20:180);
set(gca,'XtickLabel',[{'0'} labels(2:2:end)]) %reverse so it looks like Nikla's data

cat = [ 0 10 20 30 40 50 60 70 80 90  100 110 120 130 140 150 160 170;
 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180]; 
means = NaN(2,size(cat,2));
stds = NaN(2,size(cat,2));
nmp = NaN(2,size(cat,2));
for c = 1:size(cat,2);
    for r = 1:2
        ind = find(ang{r} >= cat(1,c) & ang{r} < cat(2,c));
        means(r,c) = mean(reaction_time{r}(ind));
        stds(r,c) = std(reaction_time{r}(ind));
        nmp(r,c) = length(ind);
    end
end
figure
hold on
errorbar(means(1,:),stds(1,:)./sqrt(nmp(1,:)))
errorbar(means(2,:),stds(2,:)./sqrt(nmp(2,:)),'r')
set(gca,'Xtick',1:18)
set(gca,'XtickLabel',labels) %reverse so it looks like Nikla's data
xlabel('Change in Saccade Angle')
ylabel('Time to Fixation (ms)')
legend('Predictable','Random')