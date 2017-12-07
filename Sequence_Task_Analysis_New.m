% Sequence Task analysis-All types
% written by Seth Konig January 25,2015

%%
%%---[1] Extract Eye Data from Cortex Files and Detect Fixations and Saccades with Cluster Fix---%%
% for forward bias data sets only. G did not have familirization block but
% should not effect random sequence results. The 1st 10 trials were always
% ignored in all cases

%---Data for 1 predictable and 1 random sequence---%
% sequence_cortex_files = {'PW140711.2','PW140714.2','PW140715.2','PW140716.2'};
% calibration_cortex_files = {'PW140711.1','PW140714.1','PW140715.1','PW140716.1'};
% item_letter = {'RG','RH','RI','RJ'};
% number_calibration_points = [25 25 25 25];

% sequence_cortex_files    = {'TT140922.1','TT140923.1','TT140924.1','TT140929.1'};%,'TT141002.1'};
% calibration_cortex_files = {'TT140922.2','TT140923.2','TT140924.2','TT140929.2'};%,'TT141002.2'};
% item_letter = {'RG','RH','RI','RK'};%'RL'};%purposefully removing RL
% number_calibration_points = [26 26 26 26];% 26];

% sequence_cortex_files    = {'TO150402.2','TO150403.3','TO150406.3','TO150407.2','TO150408.2','TO150409.2'};
% calibration_cortex_files = {'TO150402.1','TO150403.2','TO150406.2','TO150407.1','TO150408.1','TO150409.1'};
% item_letter = {'RG','RH','RI','RJ','RK','RL'};
% number_calibration_points = [25 25 25 25 25 25];

% sequence_cortex_files = {'RR141113.2','RR141114.2','RR141117.2','RR141118.2','RR141119.2','RR141120.2'};
% calibration_cortex_files = {'RR141113.1','RR141114.1','RR141117.1','RR141118.1','RR141119.1','RR141120.1'};
% item_letter = {'RG','RH','RJ','RI','RK','RL'};
%  number_calibration_points = 24*ones(1,6);

%Timmy attempt 2 since not enough data during attempt 1
% sequence_cortex_files = {'TT150519.2','TT150520.2','TT150521.2','TT150522.2',...
%                          'TT150526.2','TT150527.2','TT150528.2','TT150601.2',...
%                          'TT150603.2','TT150604.2','TT150605.2'};
% calibration_cortex_files =  {'TT150519.1','TT150520.1','TT150521.1','TT150522.1',...
%                              'TT150526.1','TT150527.1','TT150528.1','TT150601.1',...
%                              'TT150603.1','TT150604.1','TT150605.1'};
% item_letter = {'SeqR51','SeqR52','SeqR53','SeqR54',...
%                'SeqR55','SeqR56','SeqR57','SeqR58',...
%                'SeqR60','SeqR61','SeqR62'};
% number_calibration_points = 25*ones(1,length(item_letter));

%Red attempt 2 since not particularly good at fixating before
% sequence_cortex_files =    {'RR150522.2','RR150526.2','RR150528.2',...
%                             'RR150529.2','RR150601.2','RR150602.2','RR150604.3',...
%                             'RR150605.2','RR150608.2'};
% calibration_cortex_files = {'RR150522.1','RR150526.1','RR150528.1',...
%                             'RR150529.1','RR150601.1','RR150602.1','RR150604.2',...
%                             'RR150605.1','RR150608.1'};
% item_letter =              {'SeqR51','SeqR52','SeqR54',...
%                             'SeqR55','SeqR56','SeqR57','SeqR59',...
%                             'SeqR60','SeqR61'};
                        

%---Data for 2 predictable sequences---%
% sequence_cortex_files = {'PW140723.2','PW140813.4','PW140814.4','PW140815.3'};
% calibration_cortex_files = {'PW140723.1','PW140813.1','PW140814.1','PW140815.1'};
% item_letter = {'RS','RQ','RT','RU'};
% number_calibration_points = 25*ones(1,length(sequence_cortex_files));

% sequence_cortex_files = {'TO150410.3','TO150413.2','TO150414.2','TO150415.2','TO150416.2',...
%                          'TO150417.2','TO150420.2'};
% calibration_cortex_files = {'TO150410.2','TO150413.1','TO150414.1','TO150415.1','TO150416.1',...
%                             'TO150417.1','TO150420.1'};
% item_letter = {'RM','RN','RO','RP','RQ','RR','RS'};
% number_calibration_points = 25*ones(1,length(sequence_cortex_files));

% sequence_cortex_files = {'TT141007.1','TT141008.1','TT141013.1','TT141014.1',...
%     'TT141015.1','TT141016.1','TT141017.1'};
% calibration_cortex_files = {'TT141007.2','TT141008.2','TT141013.2','TT141014.2',...
%     'TT141015.2','TT141016.2','TT141017.2'};
% item_letter = {'RM','RN','RO','RP','RQ','RR','RS'};
% number_calibration_points = 25*ones(1,length(sequence_cortex_files));


% sequence_cortex_files = {'RR150304.2','RR150305.2','RR150306.2','RR150309.2',...
%     'RR150310.2','RR150311.2','RR150313.2'};
% calibration_cortex_files =  {'RR150304.1','RR150305.1','RR150306.1','RR150309.1',... 
%     'RR150310.1','RR150311.1','RR150313.1'};
% item_letter = {'RT','RU','RV','RW','RX','RY','RZ'};
% number_calibration_points = 24*ones(1,length(sequence_cortex_files));

% % TT reacclimation sessions
% sequence_cortex_files = {'TT150324.2','TT150325.2','TT150326.2','TT150327.2',...
%     'TT150330.2'};
% calibration_cortex_files =  {'TT150324.1','TT150325.1','TT150326.1','TT150327.1',...
%     'TT150330.1'};
% item_letter = {'RT','RU','RV','RW','RX'};
% number_calibration_points = 25*ones(1,length(sequence_cortex_files));

%TO reacclimation sessions prior to recording
% sequence_cortex_files = {'TO151105.3','TO151112.3','TO151113.3',...
%     'TO151116.4','TO151117.3','TO151119.3','TO151120.3','TO151123.3'};
% calibration_cortex_files =  {'TO151117.1','TO151117.1','TO151117.1',...
%     'TO151117.1','TO151117.1','TO151119.1','TO151120.1','TO151123.1'};
% item_letter = {'RT','RS','RU','RV','RW','RX','RY','RZ'};


%--PW Post-lesion reacclimation--%

% sequence_cortex_files = {'PW160106.2','PW160107.2','PW160108.2','PW160111.2','PW160113.2',...
%     'PW160121.2','PW160122.2'};
% calibration_cortex_files =  {'PW160106.1','PW160107.1','PW160108.1','PW160111.1','PW160113.1',...
%     'PW160121.1','PW160122.1'};
% item_letter = {'Seq07','Seq08','Seq09','Seq10','Seq12',...
%     'Seq17','Seq18'};

%--RR post-lesion reacclimation---%
% sequence_cortex_files = {'RR160119.2','RR160120.2','RR160121.2','RR160122.2',...
%                          'RR160125.2','RR160126.2'};
% calibration_cortex_files =  {'RR160119.1','RR160120.1','RR160121.1','RR160122.1',...
%                           'RR160125.1','RR160126.1'};
% item_letter = {'Seq06','Seq07','Seq08','Seq09','Seq10','Seq11'};

%%---Data for PW to determine if there is drop in predictive saccade rate when images are introduced---%
% calibration_cortex_files = {'PW150304.1','PW150305.1','PW150310.1','PW150312.1'};
% sequence_cortex_files =  {'PW150304.2','PW150305.2','PW150310.2','PW150312.2'};
% item_letter = {'RZ','SEQ01','SEQ02','SEQ03'};
% number_calibration_points = 25*ones(1,length(sequence_cortex_files));

%---Data for PW to determine 15 min of dms was a good distraction ---%
% calibration_cortex_files = {'PW150317.1','PW150317.1'};
% sequence_cortex_files =  {'PW150317.2','PW150317.4'};
% item_letter = {'SEQ04','SEQ04'};
% number_calibration_points = 25*ones(1,length(sequence_cortex_files));

%---Data for PW post-lesion random reaction times---%
% sequence_cortex_files = {'PW160222.2','PW160223.2','PW160224.2',...
%                 'PW160225.2','PW160229.2','PW160301.2','PW160302.2',...
%                 'PW160303.2','PW150304.2'};
% item_letter = [cell(1,7),'SeqR08','SeqR09'];
% calibration_cortex_files = [cell(1,7) {'PW160303.1','PW150304.1'}];
          

%---Data for RR post-lesion random reaction times---%
% sequence_cortex_files = {'RR160308.1','RR160309.1','RR160310.1','RR160311.1',...
%                          'RR160315.1','RR160316.1','RR160317.1','RR160318.1',...
%                          'RR160321.1'};
% item_letter = {'SeqR63','SeqR64','SeqR65','SeqR67',...
%                'SeqR68','SeqR69','SeqR70','SeqR71',...
%                'SeqR72'};
% calibration_cortex_files = {'RR160308.2','RR160309.2','RR160310.2','RR160311.2',...
%                             'RR160315.2','RR160316.2','RR160317.2','RR160318.2',...
%                             'RR160321.2'};

%---Data for Manfred pre-lesion random reaction times---%
% sequence_cortex_files = {'MF161104.2','MF161107.2','MF161109.3','MF161110.3','MF161117.2',...
%     'MF161121.2','MF161122.2'};
% calibration_cortex_files = {'MF161104.1','MF161107.1','MF161109.1','MF161110.1','MF161117.1',...
%     'MF161121.1','MF161122.1'};
% item_letter = {'RG','RH','RJ','RK','SeqR51',...
%     'SeqR52','SeqR53'};
%skipping RI since wasn't forward biased probably didn't load proper
%condition file; nonforward bias would produce slower RTs
%SeqR51.itm was default RL itm file not properly loaded 
          

%---Data for 2 random sequences---%
% sequence_cortex_files = {'PW140701.4','PW140702.2','PW140702.4','PW140703.2','PW140703.4'};
% calibration_cortex_files = {{'PW140701.1','PW140701.3'},'PW140702.1','PW140702.1',...
%     'PW140703.1','PW140703.1','PW140707.1'};
% item_letter = {'RA','RD','RD','RB','RB','RC'};
% number_calibration_points = [63 25 25 25 25 25];

% sequence_cortex_files = {'TO150401.2'};
% calibration_cortex_files = {'TO150401.1'};
% item_letter = {'RF'};
% number_calibration_points = [25];

% sequence_cortex_files = {'TT140915.1','TT140910.1','TT140911.2','TT140912.1','TT140915.1'};
% calibration_cortex_files = {'TT140915.1','TT140910.2','TT140911.1','TT140912.2','TT140915.2'};
% item_letter = {'RA','RB','RC','RD','RE'};
% number_calibration_points = [26 26 26 26 26 26];

%Tobii post-lesion re-acclimation--%
% sequence_cortex_files = {'TO170214.2','TO170215.2','TO170217.2','TO170221.2'};
% calibration_cortex_files = {'TO170214.1','TO170215.1','TO170217.1','TO170221.1'};
% item_letter = {'Seq38','Seq39','Seq40','Seq41','Seq42'};

%---Tobii  post-lesion random sequences---%
sequence_cortex_files = {'TO170323.2','TO170324.2','TO170327.2','TO170328.2','TO170329.2',...
                        'TO170330.2','TO170331.2'};
calibration_cortex_files =  {'TO170323.1','TO170324.1','TO170327.1','TO170328.1','TO170329.1',...
                        'TO170330.1','TO170331.1'};
item_letter = {'SeqR11','SeqR12','SeqR13','SeqR14','SeqR15','SeqR16','SeqR17'};

number_calibration_points = 25*ones(1,length(item_letter));
% for file = 1:length(sequence_cortex_files);
%     getSequenceData(sequence_cortex_files{file},calibration_cortex_files{file},...
%         item_letter{file},number_calibration_points(file));
% %     close all
% end
emailme('Done processing sequence data')
%%
%%---[2] Combine Behavior Across several different item sets---%%
eyedat_dir = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\'];
fixwin = 5; %dva of fixation window
clr = ['bycy'];

fixation_files = sequence_cortex_files;

condition = cell(1,length(fixation_files));%condition number for each trial
item_locs = cell(1,length(fixation_files));%item locations for each trial
fixation_numbers = cell(1,length(fixation_files));%fixation number associated with each item
t2f = cell(1,length(fixation_files));%time to fixation according to Cluster Fix
cortex_t2f = cell(1,length(fixation_files));%time to fixation according to Cortex
fixation_accuracy = cell(1,length(fixation_files));%accuracy of fixation on item
extrafixations = cell(1,length(fixation_files));%number of extra fixations between/before items
time_to_leave = cell(1,length(fixation_files));%time it takes monkey to saccade away from item after it disappears
fixation_duration = cell(1,length(fixation_files));%fixation duration for fixation(s) on item
cortex_predicted = cell(1,length(fixation_files));%whether cortex called it a predictive saccade
cortex_refixated = cell(1,length(fixation_files));%broke then refixated item succesffuly
which_sequence = cell(1,length(fixation_files)); %whether trial was random or predictable sequence

percentage = zeros(1,length(fixation_files)); %how well they completed the task
percentage_75 =  zeros(1,length(fixation_files)); %how well they completed the task
percentage_completed = zeros(1,length(fixation_files)); %how well they completed the task
for file = 1:length(fixation_files);
    load([eyedat_dir fixation_files{file}(1:8) '_' fixation_files{file}(end) '-fixation.mat'])
    
    if per(end).cnd == 425 && length(per) < 424 %then skipps error trials
            percentage_completed(file) = round(100*(all_trials(1,end)-1)/424);
            percentage(file) = round(100*(sum(all_trials(2,:))/(all_trials(1,end)-1)));
    else
        percentage_completed(file) = round(100*(all_trials(1,end)-1)/424);
        percentage(file) = round(100*(all_trials(1,end)-1)/size(all_trials,2));
    end
    
    condition{file} = NaN(length(fixationstats),1); %which cross hair was presented
    item_locs{file} = cell(1,length(fixationstats)); %store locations as well
    fixation_numbers{file} = NaN(length(fixationstats),4);%fixation number associated with each item
    t2f{file} = NaN(length(fixationstats),4);%time to fixation according to Cluster Fix
    cortex_t2f{file} = NaN(length(fixationstats),4);%time to fixation according to Cortex
    fixation_accuracy{file} = NaN(length(fixationstats),4);%accuracy of fixation on item
    extrafixations{file} = NaN(length(fixationstats),4);%number of extra fixations between/before items
    time_to_leave{file} = NaN(length(fixationstats),4);%time it takes monkey to saccade away from item after it disappears
    fixation_duration{file} = NaN(length(fixationstats),4);%fixation duration for fixation(s) on item
    cortex_predicted{file} = NaN(length(fixationstats),4);%whether cortex called it a predictive saccade
    cortex_refixated{file} = NaN(length(fixationstats),4);%broke then refixated item succesffuly
    which_sequence{file} = NaN(length(fixationstats),1); %whether trial was random or predictable sequence
    
%     figure 
     
    for t = 1:length(fixationstats)
        trialcnd = per(t).cnd;
        if length(find(trialcnd == all_trials(1,:))) == 1 %if only saw this sequence 1 time
            
            trialdata = analyze_sequence_trial(fixationstats{t},double(item_locations{1,trialcnd-1}),fixwin,...
                per(t).allval,per(t).alltim,false);
            
            condition{file}(t) = per(t).cnd;
            item_locs{file}{t} = item_locations{1,trialcnd-1};
            fixation_numbers{file}(t,:) = trialdata.fixationnums;
            t2f{file}(t,:) = trialdata.t2f;
            cortex_t2f{file}(t,:) = trialdata.cortext2f;
            fixation_accuracy{file}(t,:) = trialdata.accuracy;
            fixation_duration{file}(t,:) = trialdata.fixation_duration;
            cortex_predicted{file}(t,:) = trialdata.cortexpredict;
            cortex_refixated{file}(t,:) = trialdata.cortexbreak;
            which_sequence{file}(t) = item_locations{2,t}; %1 is predictable, 2 is random
            
%             if  item_locations{2,trialcnd-1} == 1
%                 subplot(1,2,1)
%                 hold on
%                 plot(fixationstats{t}.XY(1,:),fixationstats{t}.XY(2,:));
%                 hold off
%             else
%                 subplot(1,2,2)
%                 hold on
%                 plot(fixationstats{t}.XY(1,:),fixationstats{t}.XY(2,:));
%                 hold off
%             end
        end
    end
    if nanmean(nanmean(t2f{file}-cortex_t2f{file})) > 25 
        disp(['Check Condition File and Calibration for file ' ...
        fixation_files{file} ' . Cortex Times do no align with Cluster Fix'])
    end
%     subplot(1,2,1)
%     axis equal
%     xlim([-13 13])
%     ylim([-10 10])
%     subplot(1,2,2)
%     axis equal
%     xlim([-13 13])
%     ylim([-10 10])
    
end
%%
for file = 1:length(fixation_files);
    means = NaN(2,4);
    sems = NaN(2,4);
    
    means(1,:) = nanmean(t2f{file}(which_sequence{file}==1,:));
    means(2,:) = nanmean(t2f{file}(which_sequence{file}==2,:));
    sems(1,:) = nanstd(t2f{file}(which_sequence{file}==1,:))...
        ./sqrt(sum(~isnan(t2f{file}(which_sequence{file}==1,:))));
    sems(2,:) = nanstd(t2f{file}(which_sequence{file}==1,:))...
        ./sqrt(sum(~isnan(t2f{file}(which_sequence{file}==1,:))));
    
    figure
    errorbar(means',sems')
    ylabel('Time to Fixations/RT (ms)')
    xlabel('Item number')
    legend('Predictable','Random')
    
    
    
end
%%
for file = 1:length(fixation_files)-1;
    figure
    hold on
    for trial = 21:size(fixation_accuracy{file});
        if which_sequence{file}(trial) == 1
            colour = 'blue';
        else
            colour = 'red';
        end
        plot(trial,nanmean(fixation_accuracy{file}(trial,:)),['.' colour])
    end
    hold off
end


%%
all_extrafixations = [];
all_fixation_duration = [];
all_rts = [];
all_cortex_rts = [];
all_random_predict = [];
all_fixation_accuracy = [];
all_cortex_predicted = [];
for file = 1:length(fixation_files);
    all_fixation_duration = [all_fixation_duration; fixation_duration{file}(21:end,:)];
    all_extrafixations = [all_extrafixations; extrafixations{file}(21:end,:)];
    all_rts = [all_rts; t2f{file}(21:end,:)];
    all_cortex_rts = [all_cortex_rts; cortex_t2f{file}(21:end,:)];
    all_random_predict = [all_random_predict; which_sequence{file}(21:end,:)];
    all_fixation_accuracy = [all_fixation_accuracy; fixation_accuracy{file}(21:end,:)];
    all_cortex_predicted = [all_cortex_predicted; cortex_predicted{file}(21:end,:)];
end
%%
cumulative_rts = sum(all_rts(:,2:4),2);
cumulative_predictable_rts = cumulative_rts(all_random_predict == 1);
cumulative_random_rts = cumulative_rts(all_random_predict == 2);
%%
cumulative_rts = NaN(2,length(fixation_files));
for file = 1:length(fixation_files);
    cumulative_rts(1,file) = nanmean(sum(cortex_t2f{file}(which_sequence{file} == 1,2:4),2));
    cumulative_rts(2,file) = nanmean(sum(cortex_t2f{file}(which_sequence{file} == 2,2:4),2));
end
%%
random_rts = all_rts(all_random_predict == 2,:);
predict_rts = all_rts(all_random_predict == 1,:);
%%
all_random_rts = [];
all_cortex_random_rts = [];
all_cortex_predict_rts = []; 
all_predict_rts = []; 
rt_means = NaN(2,length(fixation_files));
for file = 1:length(fixation_files)
    predict_rts = t2f{file}(which_sequence{file} == 1,2:4); 
    predict_rts = predict_rts(11:end,:);
    %random_rts = t2f{file}(which_sequence{file} == 2,2:4);
    random_rts = t2f{file}(which_sequence{file} == 2,1:4);
    random_rts = random_rts(11:end,:);
    all_random_rts = [all_random_rts; random_rts(11:end,:)];
    all_predict_rts = [all_predict_rts; predict_rts];
    
    cortex_predict_rts = cortex_t2f{file}(which_sequence{file} == 1,2:4);
    cortex_predict_rts = cortex_predict_rts(11:end,:);
    %random_rts = t2f{file}(which_sequence{file} == 2,2:4);
    cortex_random_rts = cortex_t2f{file}(which_sequence{file} == 2,1:4);
    cortex_random_rts = cortex_random_rts(11:end,:);
    all_cortex_random_rts = [all_cortex_random_rts; cortex_random_rts(11:end,:)];
    all_cortex_predict_rts = [all_cortex_predict_rts; cortex_predict_rts];

    rt_means(1,file) = 100*sum(predict_rts(1:end) < 135)./sum(~isnan(predict_rts(1:end)));
    rt_means(2,file) = 100*sum(random_rts(1:end) < 135)./sum(~isnan(random_rts(1:end)));

end
%% Analysis by block (setting as every 10 trials of a given sequence)
block_prop1 = NaN(length(fixation_files),20);
block_prop2 = NaN(length(fixation_files),20);
for file = 1:length(fixation_files)
    seq1_rts = t2f{file}(which_sequence{file} == 1,2:4); 
    seq2_rts = t2f{file}(which_sequence{file} == 2,2:4); 
    
    numblock1 = floor(size(seq1_rts)/10);
    numblock2 = floor(size(seq2_rts)/10);
    
    for n1 = 1:numblock1
        temp = seq1_rts(10*(n1-1)+1:n1*10,:); 
        temp = temp(1:end); 
        block_prop1(file,n1) = sum(temp < 135)./sum(~isnan(temp)); 
    end
    
    for n2 = 1:numblock2
        temp = seq2_rts(10*(n2-1)+1:n2*10,:); 
        temp = temp(1:end); 
        block_prop2(file,n2) = sum(temp < 135)./sum(~isnan(temp)); 
    end
end
bp = (block_prop1+block_prop2)/2; 
plot(bp')
%%
ar = all_random_rts(:,2:4);
prctile_5 = prctile(ar(1:end),5);
n = sum(~isnan(ar(1:end)));
figure
hold on
hist(ar(1:end),150)
xlim([-250 450])
set(gca,'Xtick',-250:50:450)
yl = ylim;
plot([prctile_5 prctile_5],[0 yl(2)],'k--')
plot([nanmean(ar(1:end)) nanmean(ar(1:end))],[0 yl(2)],'r--')
hold on
box off
title(['Histogram of Reaction times for Items 2-4, n = ' num2str(n) ' saccades'])
xlabel('Reaction Time (ms)')
ylabel('Count')
%%

[~,pvals(1)] = ttest2(all_random_rts(:,2),all_random_rts(:,3));
[~,pvals(2)] = ttest2(all_random_rts(:,3),all_random_rts(:,4));
[~,pvals(3)] = ttest2(all_random_rts(:,2),all_random_rts(:,4));

means = nanmean(all_random_rts);
sems = nanstd(all_random_rts)./sqrt(sum(~isnan(all_random_rts))); 
figure
hold on
bar(means)
errorb(means,sems)
if pvals(1) < 0.05
    plot(2.5,max(means)+1.5*max(sems),'*k')
end
if pvals(2) < 0.05
    plot(3.5,max(means)+1.5*max(sems),'*k')
end
if pvals(3) < 0.05
    plot(3,max(means)+5.5*max(sems),'*r')
end
text(1,means(1)+2.5*sems(1),num2str(means(1)),'HorizontalAlignment','center')
text(2,means(2)+2.5*sems(2),num2str(means(2)),'HorizontalAlignment','center')
text(3,means(3)+2.5*sems(3),num2str(means(3)),'HorizontalAlignment','center')
text(4,means(4)+2.5*sems(4),num2str(means(4)),'HorizontalAlignment','center')
hold off
set(gca,'Xtick',[1:4])
xlabel('Item #')
ylabel('Reaction Times (ms)')
title('Reaction Times for Items 1-4 for Random Sequences')

%%
%plot data by session
min_rt = prctile_5;
means = zeros(1,length(t2f));
stds = zeros(1,length(t2f));
nt = zeros(1,length(t2f));
for sess = 1:length(t2f);
   t = t2f{sess}(:,2:4);
   t = t(1:end);
   means(sess) = nanmean(t);
   stds(sess) = nanstd(t);  
   nt(sess) = sqrt(size(t2f{sess},1));
end

figure
hold on
bar(means)
errorb(means,stds./nt)
xlabel('session #')
ylabel('Average RT (ms)')

% min_rt = 156;
pct = zeros(1,length(t2f));
for sess = 1:length(t2f);
    pct(sess) = sum(sum(t2f{sess}(:,2:4) < min_rt))/sum(sum(~isnan(t2f{sess}(:,2:4))));
end
xlabel('session #')
ylabel('Mean RT')

figure
bar(100*pct)
xlabel('session #')
ylabel('% of RTs < minimum rt')

%%
all_cortex_rts = [];
for sess = 1:length(cortex_t2f)
    all_cortex_rts = [all_cortex_rts; cortex_t2f{sess}];
end
