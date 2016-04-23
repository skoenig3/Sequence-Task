% Predictable 2 Sequence Task Analysis
% Written in June 2014 by Seth König adapted for 2 predictiable sequences on
% July 17, 2014
% Essentially the same as the "Sequence_Task_Analysis.m" except adatped for
% analysis for just 2 Predictable 1 sequences with many more trials

% [1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix
% [2] Combine Behavior Across several different item sets

%%
%%---[1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix---%%
% for forward bias data sets only. G did not have familirization block but
% should not effect Predictable 2 sequence results. The 1st 10 trials were always
% ignored in all cases
sequence_cortex_files = {'TT141007.1','TT141008.1','TT141013.1','TT141014.1',...
    'TT141015.1','TT141016.1','TT141017.1'};
calibration_cortex_files = {'TT141007.2','TT141008.2','TT141013.2','TT141014.2',...
    'TT141015.2','TT141016.2','TT141017.2'};

item_letter = ['MNOPQRS'];
number_calibration_points = 25*ones(1,length(sequence_cortex_files));

for file = 1:length(sequence_cortex_files);
    getSequenceData(sequence_cortex_files{file},calibration_cortex_files{file},...
        item_letter(file),number_calibration_points(file));
    %close all
end
%%
%%---[2] Combine Behavior Across several different item sets---%%
eyedat_dir = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\'];
eye_tracking_error = 5; %dva of fixation window
clr = ['bycy'];
fixation_files = sequence_cortex_files(4:7);%TT's predictive ones

IIT = cell(1,length(fixation_files));
extra_fixationpdf =  cell(1,length(fixation_files));
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
    IIT{file} = NaN(length(fixationstats),size(per(1).test,1)-1);
    
    for cross =1:4
        extrafixationpdf{file}{1,cross} = zeros(19,25);
        extrafixationpdf{file}{2,cross} = zeros(19,25);
    end
    for trial = 1:length(fixationstats)
        if  length(find(per(trial).cnd == all_trials(1,:))) == 1%otherwise saw sequence multiple times
            
            for c = 1:size(per(1).test,1)-1;
                IIT{file}(trial,c) = per(trial).test(c+1,1)-per(trial).test(c,2); 
            end
            
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
            
%                     figure
%                     hold on
%                     plot(xy(1,:),xy(2,:),'g');
%             %         for f = 1:length(fixationtimes);
%             %             plot(fixations(1,f),fixations(2,f),'k*','markersize',5);
%             %         end
%                     for c = 1:size(crosses,2);
%                         plot(crosses(1,c),crosses(2,c),'kx','markersize',12)
%                     end
%                     xlim([-17 17])
%                     ylim([-14 14])
            
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
                        if thefix ~= 1
                            extrafixations{file}(trial,c) = length(potential_fix)-1;
                            extra_fix = potential_fix;
                            extra_fix(thefix) = [];
                            for e = 1:length(extra_fix);
                                fixx = round(fixations(1,extra_fix(e)))+13;
                                fixy = round(fixations(2,extra_fix(e)))+10;
                                %remove fixations way outside of
                                %crosshair display locations
                                fixy(fixx < 1) = [];
                                fixx(fixx < 1) = [];
                                fixy(fixx > 25) = [];
                                fixx(fixx > 25) = [];
                                fixx(fixy < 1) = [];
                                fixy(fixy < 1) = [];
                                fixx(fixy > 19) = [];
                                fixy(fixy > 19) = [];
                                if ~isempty(fixx)
                                    if crosshair_locations{2,per(trial).cnd-1} == 1
                                        extrafixationpdf{file}{1,c}(fixy,fixx) = extrafixationpdf{file}{1,c}(fixy,fixx)+1;
                                    else
                                        extrafixationpdf{file}{2,c}(fixy,fixx) = extrafixationpdf{file}{2,c}(fixy,fixx)+1;
                                    end
                                end
                            end
                        end
                        fixstart = fixationtimes(1,potential_fix(thefix));
                        valid_fixations(1:potential_fix(thefix)) = 0;
                        potential_fix = potential_fix(thefix);
                        dist = dist(thefix);
                        
                    end
                    
                    if crosshair_locations{2,per(trial).cnd-1} == 1
                        seq_condition{file}(trial) = 1; %predictable 1
                    else
                        seq_condition{file}(trial) = 2;% Predictable 2
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
%                                     plot(xy(1,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),...
%                                         xy(2,fixationtimes(1,potential_fix):fixationtimes(2,potential_fix)),'c');
%                                     plot(xy(1,fixstart),xy(2,fixstart),'m*','markersize',10)
                end
            end
%                     pause(0.5)
%                     close
        end
    end
end
%%
all_reaction_time = cell(2,212); %how long after a stimulus disappears do they saccade
all_time_to_fixation = cell(2,212); %how fast did they get to it the first time around
all_fixation_accuracy = cell(2,212); %how far off
all_fixation_duration = cell(2,212); %fixation duration
all_fixation_on_test = cell(2,212); %how long the stimulus was on while fixating
all_extrafixations = cell(2,212); %how many noncross hair fixations do they make
for file = 1:length(condition);
    for cnd = 1:4
        trials = find(condition{file} == cnd);
        for repitition = 1:length(trials)
            if max(condition{file}) == 2
                if condition{file}(trials(repitition)) <= 1;
                    rand_row = 1;
                else
                    rand_row = 2;
                end
            else
                if condition{file}(trials(repitition)) <= 2;
                    rand_row = 1;
                else
                    rand_row = 2;
                end
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

if 1%length(fixation_files) > 1;
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
    for lt = 1:150
        for cnd = 1:2 %row 1 is control, row 2 is Predictable 2
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
        errorbar(mean_fixation_accuracy{1}(cross,:),std_fixation_accuracy{1}(cross,:)./sqrt(num_fixation_accuracy{1}(1,:)))
        errorbar(mean_fixation_accuracy{2}(cross,:),std_fixation_accuracy{2}(cross,:)./sqrt(num_fixation_accuracy{2}(1,:)),'r')
        hold off
        xlabel('Trial #')
        ylabel('Fixation Accuracy (dva)')
        title(['Cross ' num2str(cross)])
        legend('Predictable 1','Predictable 2')
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
        legend('Predictable 1','Predictable 2')
        ylim([0 2])
    end
    subtitle('# of Extra Fixations')
    
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
        legend('Predictable 1','Predictable 2')
    end
    subtitle('Time to Fixation')
end
%%
frt = cell(1,2);
for cnd = 1:2
    for trial = 11:length(all_time_to_fixation);
        frt{cnd} = [frt{cnd};all_time_to_fixation{cnd,trial}];
    end
end
means = NaN(2,4);
stds = NaN(2,4);
for cnd = 1:2
    for cross = 1:4
        means(cnd,cross) = nanmean(frt{cnd}(:,cross));
        stds(cnd,cross) = nanstd(frt{cnd}(:,cross));
    end
end
figure
errorb(means',stds')
ylabel('Time to Fixation (ms)')
xlabel('Cross #')
title('Mean Time to Fixation')
legend('Predictable 1','Predictable 2')

figure
hold on
bar(1,nanmean(sum(frt{1}(:,2:end),2)))
bar(2,nanmean(sum(frt{2}(:,2:end),2)),'r')
errorb(1,nanmean(sum(frt{1}(:,2:end),2)),nanstd(sum(frt{1}(:,2:end),2)))
errorb(2,nanmean(sum(frt{2}(:,2:end),2)),nanstd(sum(frt{2}(:,2:end),2)))
ylabel('Cummulative Time to Fixation (ms)')
set(gca,'Xtick',[1:2])
set(gca,'XtickLabel',{'Predictable 1','Predictable 2'})
ylim([450 700])

%%
%plot individual day data
plots = {1,2,4,5,[3 6]};
for file = 1:length(time_to_fixation);
    figure
    trials = find(seq_condition{file} == 1);
    rtrials = find(seq_condition{file} == 2);
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
        legend('Predictable 1','Predictable 2');
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
%%
% how do reaction times during Predictable 2 sequences change with trial or time
% into the task? Can't use a correlation analysis because reaction times by
% trial are not necessarily independent
t2f = [];
figure
hold on
for trial = 1:length(all_time_to_fixation);
    points = all_time_to_fixation{2,trial}(:,2:end);
    nmp = sum(sum(~isnan(points)));
    if nmp > 6 %have at least 2 trials worth of data, 3 points per trial
        mn = nanmedian(points(1:end));
        plot(trial,mn,'.')
        t2f(trial) = mn;
    end
end
ylabel('Mean Time to Fixation (ms)')
xlabel('Trail Number excluding 1st 10 trials')
title('Mean Reaction Time for Items 2-4 by Trial Across 4 sessions')
t2f(isnan(t2f)) = [];
trials = 1:length(t2f);
p = polyfit(trials,t2f,1)
plot(1:length(t2f),p(1)*(1:length(t2f))+p(2),'k')
%%
% are reaction times during Predictable 2 sequences different for different items
% in the sequence or are reaction times flat irrespective of 

[~,pval12] = ttest2(frt{2}(:,1),frt{2}(:,2))
[~,pval13] = ttest2(frt{2}(:,1),frt{2}(:,3))
[~,pval14] = ttest2(frt{2}(:,1),frt{2}(:,4))
[~,pval23] = ttest2(frt{2}(:,2),frt{2}(:,3))
[~,pval24] = ttest2(frt{2}(:,2),frt{2}(:,4))
[~,pval34] = ttest2(frt{2}(:,3),frt{2}(:,4))

f =frt{2}(:,2:end);
prct5 = prctile(f(1:end),5)

figure
hist(f(1:end),50);
hold on
plot([nanmean(f(1:end)) nanmean(f(1:end))],[0 350],'k')
plot([prct5 prct5],[0 350],'r')
xlabel('Reaction time (ms)')
ylabel('Count')
%%
extrapdf = zeros(19,25);
for i = 1:numel(extrafixationpdf{1});
    extrapdf =  extrapdf+extrafixationpdf{1}{i};
end