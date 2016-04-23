function [fixationtimes,fixations] = remove_corrective_saccades(...
    xy,saccadetimes,saccade_amplitudes,eye_tracking_error)
% written by Seth Konig on June 26, 2014
% function removes "corrective saccades" determine as saccades (other than
% the 1st recorded saccade) that are shorter than a defined threshold.
% This threshold is set by the "saccade_amplitude" variable.

corrective_saccades = find(saccade_amplitudes <= eye_tracking_error/2);
corrective_saccades(corrective_saccades == 1) = []; %ignore these
saccadetimes(:,corrective_saccades) =[];%remove corrective saccades

saccade_indexes = [];
for sac = 1:size(saccadetimes,2)
    saccade_indexes = [saccade_indexes saccadetimes(1,sac):saccadetimes(2,sac)];
end

fixation_indexes = 1:size(xy,2);
[~ , ia, ~] = intersect(fixation_indexes,saccade_indexes);
fixation_indexes(ia) = [];

[fixationtimes,fixations] = BehavioralIndexXY(fixation_indexes,xy(1,:),xy(2,:));

    function [behaviortime, behaviormean] = BehavioralIndexXY(behavind,x,y)
        %function is the same as above but also calculates mean fixation position
        dind = diff(behavind);
        gaps =find(dind > 1);
        behaveind = zeros(length(gaps),50);
        if ~isempty(gaps)
            for gapind = 1:length(gaps)+1;
                if gapind == 1;
                    temp = behavind(1:gaps(gapind));
                elseif gapind == length(gaps)+1
                    temp = behavind(gaps(gapind-1)+1:end);
                else
                    temp = behavind(gaps(gapind-1)+1:gaps(gapind));
                end
                behaveind(gapind,1:length(temp)) = temp;
            end
        else
            behaveind =  behavind;
        end
        behaviortime = zeros(2,size(behaveind,1));
        behaviormean = zeros(2,size(behaveind,1));
        for index=1:size(behaveind,1)
            rowfixind = behaveind(index,:);
            rowfixind(rowfixind == 0) = [];
            behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
            behaviormean(:,index) = [mean(x(rowfixind));...
                mean(y(rowfixind))];
        end
    end
end