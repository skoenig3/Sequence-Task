% Sequence Task Analysis
% Written in June 2014 by Seth König
% [1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix
% [2] Combine Behavior Across several different item sets

%%
%%---[1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix---%%
sequence_cortex_files = {'TT140707.4'};
calibration_cortex_files = {'TT140707.3'};
item_number = [7];
number_calibration_points = [25];

for file = 1:length(sequence_cortex_files);
    getSequenceData(sequence_cortex_files{file},calibration_cortex_files{file},...
        item_number(file),number_calibration_points(file));
    close all
end
%%
%%---[2] Combine Behavior Across several different item sets---%%
eyedat_dir = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\'];
eye_tracking_error = 5; %dva of fixation window
clr = ['bycy'];

files = dir(eyedat_dir);
fixation_files = [];
for f = 1:length(files);
    if ~isempty(strfind(files(f).name,'-fixation.mat'))
        load([eyedat_dir files(f).name],'itemnum');
        if ~ischar(itemnum)
            if strcmpi(files(f).name(1:2),'TT')
                fixation_files = [fixation_files {files(f).name}];
            end
        end
    end
end

condition = cell(1,length(fixation_files));
reaction_time = cell(1,length(fixation_files)); %how long after a stimulus disappears do they saccade
time_to_fixation = cell(1,length(fixation_files)); %how fast did they get to it the first time around
fixation_accuracy = cell(1,length(fixation_files)); %how far off
fixation_duration = cell(1,length(fixation_files)); %fixation duration
fixation_on_test = cell(1,length(fixation_files)); %how long the stimulus was on while fixating
extrafixations = cell(1,length(fixation_files)); %how many noncross hair fixations do they make
for file = 1:length(fixation_files);
    load([eyedat_dir fixation_files{file}])
    condition{file} = NaN(length(fixationstats),1); %which cross hair was presented
    reaction_time{file} = NaN(length(fixationstats),size(per(1).test,1)-1); %how long after a stimulus disappears do they saccade
    time_to_fixation{file} = NaN(length(fixationstats),size(per(1).test,1)); %how fast did they get to it the first time around
    fixation_accuracy{file} = NaN(length(fixationstats),size(per(1).test,1)); %how far off
    fixation_duration{file} = NaN(length(fixationstats),size(per(1).test,1)); %fixation duration
    fixation_on_test{file} = NaN(length(fixationstats),size(per(1).test,1)); %how long the stimulus was on while fixating
    extrafixations{file} = zeros(length(fixationstats),size(per(1).test,1)); %how many noncross hair fixations do they make
    
    
    successful_trials = find(all_trials(2,:) == 1);
    trial = 1:length(fixationstats);
    bad_trials = [];
    for trials = 1:length(fixationstats)
        if trials == 1;
            if all_trials(2,1) == 0;
                bad_trials = [bad_trials 1];
            end
        else
            if all_trials(2,successful_trials(trials)-1) == 0; %previous trial was not a failure
                bad_trials = [bad_trials trials];
            end
        end
    end
    %everything but the 1st presentation of sequence doens't really matter here
    bad_trials(bad_trials > 1 & bad_trials<= 10) = [];
    bad_trials(bad_trials > 11 & bad_trials<= 20) = [];
    bad_trials(bad_trials > 21 & bad_trials<= 30) = [];
    bad_trials(bad_trials > 31 & bad_trials<= 40) = [];
    
    for trial = 1:length(fixationstats)
        if isempty(find(bad_trials == trial))
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
            
            crosses = double(crosshair_locations{per(trial).cnd}(:,1:size(per(trial).test,1)));
            
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
                        fixstart = fixationtimes(1,potential_fix(thefix));
                        valid_fixations(1:potential_fix(thefix)) = 0;
                        potential_fix = potential_fix(thefix);
                        dist = dist(thefix);
                        if thefix ~= 1
                            extrafixations{file}(trial,c) = thefix-1;
                        end
                    end
                    
                    condition{file}(trial) = per(trial).cnd;
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
        else
              condition{file}(trial) = per(trial).cnd;
        end
    end
end

all_reaction_time = cell(1,50); %how long after a stimulus disappears do they saccade
all_time_to_fixation = cell(1,50); %how fast did they get to it the first time around
all_fixation_accuracy = cell(1,50); %how far off
all_fixation_duration = cell(1,50); %fixation duration
all_fixation_on_test = cell(1,50); %how long the stimulus was on while fixating
all_extrafixations = cell(1,50); %how many noncross hair fixations do they make
for file = 1:length(condition);
    for cnd = 1:4
        trials = find(condition{file} == cnd);
        for repitition = 1:length(trials)
            all_reaction_time{repitition} = [all_reaction_time{repitition}; reaction_time{file}(trials(repitition),:)];
            all_time_to_fixation{repitition} = [all_time_to_fixation{repitition}; time_to_fixation{file}(trials(repitition),:)];
            all_fixation_accuracy{repitition} = [all_fixation_accuracy{repitition}; fixation_accuracy{file}(trials(repitition),:)];
            all_fixation_duration{repitition} = [all_fixation_duration{repitition}; fixation_duration{file}(trials(repitition),:)];
            all_fixation_on_test{repitition} = [all_fixation_on_test{repitition}; fixation_on_test{file}(trials(repitition),:)];
            all_extrafixations{repitition} = [all_extrafixations{repitition}; extrafixations{file}(trials(repitition),:)];
        end
    end
end

too_few = cellfun(@numel,all_reaction_time);
too_few = find(too_few < 10);
all_reaction_time(too_few) =[];
all_time_to_fixation(too_few)=[];
all_fixation_accuracy(too_few)=[];
all_fixation_duration(too_few)=[];
all_fixation_on_test(too_few)=[];
all_extrafixations(too_few)=[];

mean_reaction_time = NaN(3,10);
mean_time_to_fixation = NaN(3,10);
mean_fixation_accuracy = NaN(4,10);
mean_fixation_duration = NaN(4,10);
mean_fixation_on_test = NaN(4,10);
mean_extrafixations = NaN(4,10);
std_reaction_time = NaN(3,10);
std_time_to_fixation = NaN(4,10);
std_fixation_accuracy = NaN(4,10);
std_fixation_duration = NaN(4,10);
std_fixation_on_test = NaN(4,10);
std_extrafixations = NaN(4,10);
num_reaction_time = NaN(3,10);
num_time_to_fixation = NaN(4,10);
num_fixation_accuracy = NaN(4,10);
num_fixation_duration = NaN(4,10);
num_fixation_on_test = NaN(4,10);
num_extrafixations = NaN(4,10);
for lt = 1:length(all_reaction_time);
    for cross = 1:4
        mean_time_to_fixation(cross,lt) = nanmean(all_time_to_fixation{lt}(:,cross));
        std_time_to_fixation(cross,lt) = nanstd(all_time_to_fixation{lt}(:,cross));
        num_time_to_fixation(cross,lt) = sum(~isnan(all_time_to_fixation{lt}(:,cross)));
        mean_fixation_accuracy(cross,lt) = nanmean(all_fixation_accuracy{lt}(:,cross));
        std_fixation_accuracy(cross,lt) = nanstd(all_fixation_accuracy{lt}(:,cross));
        num_fixation_accuracy(cross,lt) = sum(~isnan(all_fixation_accuracy{lt}(:,cross)));
        mean_fixation_on_test(cross,lt) = nanmean(all_fixation_on_test{lt}(:,cross));
        std_fixation_on_test(cross,lt) = nanstd(all_fixation_on_test{lt}(:,cross));
        num_fixation_on_test(cross,lt) = sum(~isnan(all_fixation_on_test{lt}(:,cross)));
        mean_extrafixations(cross,lt) = nanmean(all_extrafixations{lt}(:,cross));
        std_extrafixations(cross,lt) = nanstd(all_extrafixations{lt}(:,cross));
        num_extrafixations(cross,lt) = sum(~isnan(all_extrafixations{lt}(:,cross)));
        mean_fixation_duration(cross,lt) = nanmean(all_fixation_duration{lt}(:,cross));
        std_fixation_duration(cross,lt) = nanstd(all_fixation_duration{lt}(:,cross));
        num_fixation_duration(cross,lt) = sum(~isnan(all_fixation_duration{lt}(:,cross)));
        if cross <= 3
            mean_reaction_time(cross,lt) = nanmean(all_reaction_time{lt}(:,cross));
            std_reaction_time(cross,lt) = nanstd(all_reaction_time{lt}(:,cross));
            num_reaction_time(cross,lt) = sum(~isnan(all_reaction_time{lt}(:,cross)));
        end
    end
end
%%
figure
errorbar(mean_time_to_fixation',std_time_to_fixation'./sqrt(num_time_to_fixation'))
xlabel('Trial #')
ylabel('Time to Fixation (ms)')
legend('Cross 1','Cross 2','Cross 3','Cross 4')

figure
errorbar(mean_fixation_accuracy',std_fixation_accuracy'./sqrt(num_fixation_accuracy'))
xlabel('Trial #')
ylabel('Fixation accuracy (dva)')
legend('Cross 1','Cross 2','Cross 3','Cross 4')

figure
errorbar(mean_fixation_on_test',std_fixation_on_test'./sqrt(num_fixation_on_test'))
xlabel('Trial #')
ylabel('Fixation Duration while on Test (ms)')
legend('Cross 1','Cross 2','Cross 3','Cross 4')

figure
errorbar(mean_fixation_duration',std_fixation_duration'./sqrt(num_fixation_duration'))
xlabel('Trial #')
ylabel('Fixation Duration/Time spent in Fixation Window(ms)')
legend('Cross 1','Cross 2','Cross 3','Cross 4')

figure
errorbar(mean_extrafixations',std_extrafixations'./sqrt(num_extrafixations'))
xlabel('Trial #')
ylabel('Number of Extra Fixations before fixating on cross')
legend('Cross 1','Cross 2','Cross 3','Cross 4')

figure
errorbar(mean_reaction_time',std_reaction_time'./sqrt(num_reaction_time'))
xlabel('Trial #')
ylabel('Reaction Time (ms)')
legend('Cross 1','Cross 2','Cross 3')

frt = [];
for trial = 11:length(all_time_to_fixation);
    frt = [frt;all_time_to_fixation{trial}];
end
meanfrt = nanmean(frt);
semfrt = nanstd(frt)./sqrt(size(frt,1)-1);

figure
hold on
bar(meanfrt);
errorb(meanfrt,semfrt);
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Cross 1','Cross 2','Cross 3','Cross 4'})
ylabel('Time to Fixation (ms)')


rt = [];
for trial = 1:length(all_reaction_time);
    rt = [rt;all_reaction_time{trial}];
end
meanrt = nanmean(rt);
semrt = nanstd(rt)./sqrt(size(rt,1)-1);

figure
hold on
bar(meanrt);
errorb(meanrt,semrt);
set(gca,'Xtick',1:3)
set(gca,'XtickLabel',{'Cross 1','Cross 2','Cross 3'})
ylabel('Reaction Time (ms)')
%%
%plot individual day data
plots = {1,2,4,5,[3 6]};
for file = 1:length(time_to_fixation);
    figure
    for cnd = 1:4
        trials = find(condition{file} == cnd);
        cnd_data = time_to_fixation{file}(trials,:);
        cnd_data(cnd_data > 700) = 725;%for plotting purposes need to cap
        cnd_data(cnd_data < -200) = -225; %for plotting purposes need to cap
        subplot(2,3,plots{cnd})
        plot(cnd_data);
        ylim([-250 750])
        xlabel('Trial #')
        ylabel('Time to Fixation')
        legend('Cross 1','Cross 2','Cross 3','Cross 4');
        title(['Sequence #' num2str(cnd)])
    end
    subplot(2,3,plots{end})
    hold on
    bar(nanmean(time_to_fixation{file}))
    errorb(nanmean(time_to_fixation{file}),nanstd(time_to_fixation{file}))
    hold off
    ylim([0 350])
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Cross 1','Cross 2','Cross 3','Cross 4'})
    ylabel('Time to Fixation (ms)')
    subtitle(fixation_files{file}(1:8))
end
%% The stupid test if calibration doesn't work out
dur = [];
conditions = [];
for trial = 1:length(per)
    dur = [dur; diff(per(trial).test')-520];
    conditions = [conditions per(trial).cnd];
end
    