% Written by Seth Konig January 28, 2015
% trying determine when the fixation starts compared to cortex and Cluster
% Fix. I want to know how sensitive Cluster Fix is to the onset of the
% fixation. Sequence task provides a good control task since there's
% saccades to a fixation target with a fixation window so we know where and
% reasonably well when they fixated. Cluster Fix says that fixations on
% average start 13 ms after Cortex registers a fixation within the fixation
% window for PW's data from the sequence task with 1 random and 1
% predictable sequence.

% also want to look at velcoity trace for glisades

sequence_cortex_files = {'PW140711.2','PW140714.2','PW140715.2','PW140716.2'};
fixation_files = sequence_cortex_files;

eyedat_dir = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task\Eye Data\'];
fixwin = 5; %dva of fixation window
clr = ['bycy'];

%---Important Cortex Codes---%
eye_data_start_code = 100; %when xy eye data can be aligned to event codes, usually off by ~1000 ms
item_on_off_codes = [23 25 27 29; 24 26 28 30];%row 1, item on, row 2 item off, columb by item#
fixation_code = 8;%monkey fixated an item

distance_from_fixation = cell(2,length(fixation_files));%distance over time
% from mean fixation location locked to "fixation" onset

distance_from_item = cell(2,length(fixation_files));%distance over time from 
%item location locked to "fixation" onset

velocity_locked_saccades = cell(2,length(fixation_files));%velocity over time 
% locked to fixation onset. Looking for evidence for or against glissades. Going 
% align to cortex and Cluster Fix "fixation onsets" in case there's a
% difference

prepostdata = 50;%number of samples before and after fixation onset to grab
for file = 1:length(fixation_files);
    load([eyedat_dir fixation_files{file}(1:8) '_' fixation_files{file}(end) '-fixation.mat'])
    
    distance_from_fixation{1,file} = NaN(length(fixationstats),2*prepostdata+1);%for Cluster Fix
    distance_from_fixation{2,file} = NaN(length(fixationstats),2*prepostdata+1);%for Cortex
    distance_from_item{1,file} = NaN(length(fixationstats),2*prepostdata+1);%for Cluster Fix
    distance_from_item{2,file} = NaN(length(fixationstats),2*prepostdata+1);%for Cortex
    
    for t = 1:length(fixationstats)
        if any(per(t).allval == 3)  %successful/rewarded trials
            % && crosshair_locations{2,t} > 1; %if want to seperate random
            % (> 1) vs predictable (1) sequences 
            item_locations = double(crosshair_locations{1,t});
            event_codes = per(t).allval;
            event_times = per(t).alltim;
            %---get the important eyedat from Cluster Fix's output---%
            fixations = fixationstats{t}.fixations;
            fixationtimes = fixationstats{t}.fixationtimes;
            saccadetimes = fixationstats{t}.saccadetimes;
            xy = fixationstats{t}.XY;
            velx = diff(xy(1,:));
            vely = diff(xy(2,:));
            velocity = sqrt(velx.^2+vely.^2);
            
            eye_data_start = event_times(event_codes == eye_data_start_code)-1;%-1 so index starts at 1 not 0
            
            %---determine which fixation was for which item in the sequence---%
            fixation_number_on_item = NaN(1,4); %the ordinal fixation number on that item
            fixation_accuracy = NaN(1,4);%how close was the fixation to the center of the item
            for item = 1:4
                if item == 1;
                    valid_window = [1 event_times(event_codes == item_on_off_codes(2,item))-eye_data_start];
                else
                    valid_window = [event_times(event_codes == item_on_off_codes(2,item-1))...
                        event_times(event_codes == item_on_off_codes(2,item))]-eye_data_start;
                end
                if valid_window(2) > length(xy)
                    disp('error: probably have an indexing issue')
                end
                
                %find fixations that occur within a time period between the last item
                %on (or trial start if 1st item) and that item turning off
                potential_fixes = find(fixationtimes(1,:) >= valid_window(1) & fixationtimes(1,:) < valid_window(2));
                if isempty(potential_fixes)
                    disp('No potentail fixation found for this item')
                end
                
                %find valid/potentail fixation inside the fixation window
                iswithin = NaN(1,length(potential_fixes)); %determine if fixation are within fixwin
                dist_from_item = NaN(1,length(potential_fixes)); %determine distance from fixation to item
                for f = 1:length(potential_fixes)
                    if fixations(1,potential_fixes(f)) >= item_locations(1,item)-fixwin/2 && ...
                            fixations(1,potential_fixes(f)) <= item_locations(1,item)+fixwin/2 && ...
                            fixations(2,potential_fixes(f)) >= item_locations(2,item)-fixwin/2 && ...
                            fixations(2,potential_fixes(f)) <= item_locations(2,item)+fixwin/2
                        iswithin(f) = 1;
                    end
                    dist_from_item(f) = sqrt((fixations(1,potential_fixes(f))-item_locations(1,item)).^2 ...
                        +(fixations(2,potential_fixes(f))-item_locations(2,item)).^2);
                end
                if ~any(iswithin) %if no fixation found
                    if any(dist_from_item < 3) %in case calibration issue which there shouldn't
                        %give an extra 0.5 dva
                        iswithin(dist_from_item < 3) = 1;
                    else
                        disp('No Fixations on Item found')
                        %not sure code and calibration is 100% accurate so don't want to
                        %call an error
                        continue;%go to the next item
                    end
                end
                
                %the first fixation in the fixation window is the right one
                the_fixation = find(iswithin == 1);
                fixation_number_on_item(item) = potential_fixes(the_fixation(1));%go back to original index
                fixation_accuracy(item) = dist_from_item(the_fixation(1));
            end
            
            for item = 1:4
                if item == 1;
                    relevant_code_indexes = [find(event_codes == item_on_off_codes(1,item))...
                        find(event_codes == item_on_off_codes(2,item))];
                else
                    relevant_code_indexes = [find(event_codes == item_on_off_codes(2,item-1))...
                        find(event_codes == item_on_off_codes(2,item))];
                end
                relevant_events = event_codes(relevant_code_indexes(1):relevant_code_indexes(2));
                relevant_event_times = event_times(relevant_code_indexes(1):relevant_code_indexes(2))...
                    -eye_data_start;%correct for offset from eye data start
                if item == 1
                    event_start = relevant_event_times(1);%when item turned on
                else
                    event_start = relevant_event_times(2);
                    %when item turned on, because 1st event is when last item turned off
                end
                
                event_end   = relevant_event_times(end);%when item turned off
                
                
                %---Get some information from Cortex Codes---%
                % 1) determine if monkey broke fixation and refixated the item
                % 2) determine when cortex thought the monkey fixated the item
                % 3) determine if cortex encoded a truly predictive saccade causing the
                % item to appear early or would have (i.e. between 0-50 ms early)
                cortexfixations = relevant_event_times(relevant_events == fixation_code);
                cortexfixations = cortexfixations(1);%incase there are refixations
                
                if  ~isnan(fixation_number_on_item)
                    fixationstart = fixationtimes(1,fixation_number_on_item(item));
                    
                    if cortexfixations-prepostdata > 0 && ...
                            cortexfixations+prepostdata < length(xy) ...
                            && fixationstart-prepostdata > 0 && ...
                            fixationstart+prepostdata < length(xy)
                        
                        %---get distance from mean fixation location according to cortex---%
                        cortex_fix_distance = sqrt(...
                            (xy(1,cortexfixations-prepostdata:cortexfixations+prepostdata)-...
                            fixations(1,fixation_number_on_item(item))).^2+...
                            (xy(2,cortexfixations-prepostdata:cortexfixations+prepostdata)-...
                            fixations(2,fixation_number_on_item(item))).^2);
                        
                        ClusterFix_fix_distance = sqrt(...
                            (xy(1,fixationstart-prepostdata:fixationstart+prepostdata)-...
                            fixations(1,fixation_number_on_item(item))).^2+...
                            (xy(2,fixationstart-prepostdata:fixationstart+prepostdata)-...
                            fixations(2,fixation_number_on_item(item))).^2);
                        
                        distance_from_fixation{1,file}(t,:) = ClusterFix_fix_distance;
                        distance_from_fixation{2,file}(t,:) = cortex_fix_distance;
                        
                        %---get distance from item location---%
                        cortex_item_distance = sqrt(...
                            (xy(1,cortexfixations-prepostdata:cortexfixations+prepostdata)-...
                            item_locations(1,item)).^2+...
                            (xy(2,cortexfixations-prepostdata:cortexfixations+prepostdata)-...
                            item_locations(2,item)).^2);
                        
                        ClusterFix_item_distance = sqrt(...
                            (xy(1,fixationstart-prepostdata:fixationstart+prepostdata)-...
                            item_locations(1,item)).^2+...
                            (xy(2,fixationstart-prepostdata:fixationstart+prepostdata)-...
                            item_locations(2,item)).^2);
                        
                        distance_from_item{1,file}(t,:) = ClusterFix_item_distance;
                        distance_from_item{2,file}(t,:) = cortex_item_distance;
                        
                        %---Lock velocity to Fixation Onset---%
                        %aka saccade offset 
                        velocity_locked_saccades{1,file}(t,:) = ...
                            velocity(fixationstart-prepostdata:fixationstart+prepostdata);
                       velocity_locked_saccades{2,file}(t,:) = ...
                            velocity(cortexfixations-prepostdata:cortexfixations+prepostdata);
                    end
                end
            end
        end
    end
end

all_cortex_distance = [];
all_ClusterFix_distance = [];
all_cortex_item_distance = [];
all_ClusterFix_item_distance = [];
all_cortex_velocites = [];
all_ClusterFix_velocities = []; 
for f = 1:length(fixation_files)
    all_ClusterFix_distance = [all_ClusterFix_distance; distance_from_fixation{1,file}];
    all_cortex_distance = [all_cortex_distance; distance_from_fixation{2,file}];
    
    all_ClusterFix_item_distance = [all_ClusterFix_item_distance; distance_from_item{1,file}];
    all_cortex_item_distance = [all_cortex_item_distance; distance_from_item{2,file}];
    
    all_ClusterFix_velocities = [all_ClusterFix_velocities; velocity_locked_saccades{1,file}];
    all_cortex_velocites = [all_cortex_velocites; velocity_locked_saccades{2,file}];

end

figure
hold on
plot(nanmean(all_cortex_distance),'r')
plot(nanmean(all_ClusterFix_distance))
hold off
xlabel('Time from "Fixation" (ms)')
ylabel('Distance from center of fixation (dva)')
set(gca,'Xtick',[0:5:2*(prepostdata)+1])
set(gca,'Xticklabel',num2cell([0:5:2*(prepostdata)+1]-prepostdata))
legend('Fixation Onset According to Cortex','Fixation Onset According to Cluster Fix')
xlim([0 2*(prepostdata+1)])
title('Distance from mean fixation location according to Cluster Fix')

figure
hold on
plot(nanmean(all_cortex_item_distance),'r')
plot(nanmean(all_ClusterFix_item_distance))
hold off
xlabel('Time from "Fixation" (ms)')
ylabel('Distance from item location (dva)')
set(gca,'Xtick',[0:5:2*(prepostdata)+1])
set(gca,'Xticklabel',num2cell([0:5:2*(prepostdata)+1]-prepostdata))
legend('Fixation Onset According to Cortex','Fixation Onset According to Cluster Fix')
xlim([0 2*(prepostdata+1)])
title('Distance from item location')

figure
hold on
plot(1000*nanmean(all_cortex_velocites),'r')
plot(1000*nanmean(all_ClusterFix_velocities))
xlabel('Time from "Fixation" (ms)')
ylabel('Velocity in dva/sec')
set(gca,'Xtick',[0:5:2*(prepostdata)+1])
set(gca,'Xticklabel',num2cell([0:5:2*(prepostdata)+1]-prepostdata))
legend('Fixation Onset According to Cortex','Fixation Onset According to Cluster Fix')
xlim([0 2*(prepostdata+1)])
title('Velocity over time locked to Fixation Onset')