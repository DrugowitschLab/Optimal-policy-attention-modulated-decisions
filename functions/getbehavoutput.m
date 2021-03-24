function outputstruct = getbehavoutput(dstruct)
%Given a data structure with behavioral data from multiple subjects,
%outputs data organized for plotting for krahbich et al., 2010 paper
%figures
warning off;

% all unique valuse/value differences
valunique = []; valdiffunique = [];
for s = 1:length(dstruct)
    valunique = cat(1,valunique,[dstruct(s).itemval(:,1);dstruct(s).itemval(:,2)]);
    valdiffunique = cat(1,valdiffunique,dstruct(s).itemval(:,1)-dstruct(s).itemval(:,2));
end
valunique = unique(valunique)';
valdiffunique = unique(valdiffunique)';
valdiffunique_abs = unique(abs(valdiffunique));

% bin specs
edges_item1fixadv = -1.5:0.2:1.5;
means_item1timeadv = edges_item1fixadv(1:end-1) + mean(diff(edges_item1fixadv))/2;
edges_firstfixdur = 0:0.2:1;
% edges_firstfixdur = 0.05:0.1:3.05;
means_firstfixdur = edges_firstfixdur(1:end-1) + mean(diff(edges_firstfixdur))/2;

% Initialize variables

% Variables for psychometric curve
rt_valdiff_abs = nan(length(dstruct),length(valdiffunique_abs));
switchcount_valdiff_abs = nan(length(dstruct),length(valdiffunique_abs));
switchcountRT_valdiff_abs = nan(length(dstruct),length(valdiffunique_abs));
rt_mean = nan(length(dstruct),1);

% Middle fixation variables
mid_fixitemval = nan(length(dstruct),length(valunique));
mid_fixunfixdiff = nan(length(dstruct),length(valdiffunique));
mid_absvaldiff = nan(length(dstruct),length(valdiffunique_abs));

% Choice biase variables
item1chosen_valdiff_all = nan(length(dstruct),length(valdiffunique));
item1chosen_valdiff_lastL = nan(length(dstruct),length(valdiffunique));
item1chosen_valdiff_lastR = nan(length(dstruct),length(valdiffunique));
item1chosen_timeadv_norm = nan(length(dstruct),length(means_item1timeadv));  % proportion of item1 chosen, by how much longer item1 was attended to
item1chosen_timeadv_beta = nan(length(dstruct),1);
firstchosen_firstfixdur = nan(length(dstruct),length(means_firstfixdur));  % proportion of first item chosen, by first fixation duration
firstchosen_firstfixdur_norm = nan(length(dstruct),length(means_firstfixdur));

% Fixation difference per value difference
fixdiff_valdiff = nan(length(dstruct),length(valdiffunique));


for s = 1:length(dstruct)
    % Variables that need to be extracted from trial-to-trial data
    firstchosen = nan(length(dstruct(s).trialnum),1);
    firstfixdur = nan(length(dstruct(s).trialnum),1);
    lastfixitem = nan(length(dstruct(s).trialnum),1);
    midfixdur_sub = [];
    mid_fixitemval_sub = []; mid_unfixitemval_sub = [];
    switchcount = nan(length(dstruct(s).trialnum),1);
    switchcount_rtNorm = nan(length(dstruct(s).trialnum),1);    % Number of switches per 1s of RT
    for ti = 1:length(dstruct(s).trialnum)
        % last fixated item
        lastfixitem(ti) = dstruct(s).fixitem{ti}(end);
        % first fixed item & duration
        firstfixitem = dstruct(s).fixitem{ti}(1);
        firstchosen(ti) = dstruct(s).choice(ti)==firstfixitem;
        firstfixdur(ti) = dstruct(s).fixdur{ti}(1);
        
        % middle fixation properties
        for fix_i = 2:length(dstruct(s).fixitem{ti})-1
            fixitem = dstruct(s).fixitem{ti}(fix_i);
            fixdur = dstruct(s).fixdur{ti}(fix_i);
            itemvals = dstruct(s).itemval(ti,:);
            midfixdur_sub = cat(1,midfixdur_sub,fixdur);
            mid_fixitemval_sub = cat(1,mid_fixitemval_sub,itemvals(fixitem));
            mid_unfixitemval_sub = cat(1,mid_unfixitemval_sub,itemvals(3-fixitem));
        end
        
        numSwitches = length(dstruct(s).fixitem{ti})-1;
        switchcount(ti) = numSwitches;
        switchcount_rtNorm(ti) = numSwitches/dstruct(s).rt(ti);
    end
    valdiff = dstruct(s).itemval(:,1)-dstruct(s).itemval(:,2);
    fixdiff = dstruct(s).tItem(:,1)-dstruct(s).tItem(:,2);
    firstchosen_valdiff = nan(1,length(valdiffunique_abs));
    item1chosen = double(dstruct(s).choice==1);
    
    % item1 fixation advandage
    item1fixadv = dstruct(s).tItem(:,1)-dstruct(s).tItem(:,2);
    
    % Get mean RT
    rt_mean(s) = mean(dstruct(s).rt);
    
    % 1. Get proportion of item1 item chosen per unique value difference
    % 2. Get relative fixation duration based on value difference
    item1chosen_norm = nan(size(item1chosen));
    for i = 1:length(valdiffunique)
        v_i = valdiff==valdiffunique(i);
        if any(v_i)
            item1chosen_valdiff_all(s,i) = mean(dstruct(s).choice(v_i)==1);
            item1chosen_valdiff_lastL(s,i) = mean(dstruct(s).choice(v_i & lastfixitem==1)==1);
            item1chosen_valdiff_lastR(s,i) = mean(dstruct(s).choice(v_i & lastfixitem==2)==1);
            % Normalized choices
            item1chosen_norm(v_i) = item1chosen(v_i) - item1chosen_valdiff_all(s,i);
            
            fixdiff_valdiff(s,i) = mean(fixdiff(v_i));
        end
    end
    % Get proportion of first item chosen per difficulty level
    firstchosen_norm = nan(size(firstchosen));
    for i = 1:length(valdiffunique_abs)
        v_i = abs(valdiff)==valdiffunique_abs(i);
        if any(v_i)
            firstchosen_valdiff(i) = mean(firstchosen(v_i));
            % Normalized choices
            firstchosen_norm(v_i) = firstchosen(v_i) - firstchosen_valdiff(i);
        end
    end
    % item1 chosen | item1 time advantage
    bin_i = discretize(item1fixadv,edges_item1fixadv);
    for b = 1:length(means_item1timeadv)
        item1chosen_timeadv_norm(s,b) = nanmean(item1chosen_norm(bin_i==b));
    end
    % first chosen | first fixation duration
    bin_i = discretize(firstfixdur,edges_firstfixdur);
    for b = 1:length(means_firstfixdur)
        firstchosen_firstfixdur(s,b) = nanmean(firstchosen(bin_i==b));
        firstchosen_firstfixdur_norm(s,b) = nanmean(firstchosen_norm(bin_i==b));
    end
    
    % middle fixation duration depending on values of the two items
    for vi = 1:length(valunique)
        mid_fixitemval(s,vi) = mean(midfixdur_sub(mid_fixitemval_sub==valunique(vi)));
    end
    for vi = 1:length(valdiffunique)
        mid_fixunfixdiff(s,vi) = mean(midfixdur_sub(mid_fixitemval_sub-mid_unfixitemval_sub==valdiffunique(vi)));
    end
    for vi = 1:length(valdiffunique_abs)
        mid_absvaldiff(s,vi) = mean(midfixdur_sub(abs(mid_fixitemval_sub-mid_unfixitemval_sub)==valdiffunique_abs(vi)));
        
        % Get RT vs. valdiff_abs for psychometric curve
        rt_valdiff_abs(s,vi) = mean(dstruct(s).rt(abs(valdiff)==valdiffunique_abs(vi)));
        switchcount_valdiff_abs(s,vi) = mean(switchcount(abs(valdiff)==valdiffunique_abs(vi)));
        switchcountRT_valdiff_abs(s,vi) = mean(switchcount_rtNorm(abs(valdiff)==valdiffunique_abs(vi)));
    end
    
    % Quantify fixation bias 
    % 1. using logistic regression
    % For choice vs. fixation advantage, fit a logistic regression curve
    [b,~,stats] = glmfit(item1fixadv,double(dstruct(s).choice==1),'binomial');
    % If coefficients perfectly separate the two conditions, b is not
    % finite. Ignore these.
    if any(stats.se > 50), b(2) = nan; end
    item1chosen_timeadv_beta(s) = b(2);
    
    % 2. using linear regression of plot
% %     y = item1chosen_timeadv_norm(s,:)';
% %     x = [ones(length(means_item1timeadv),1),means_item1timeadv'];
% %     b = regress(y,x);
% %     item1chosen_timeadv_beta(s) = b(2);
end

outputstruct = struct;
outputstruct.val = valunique;
outputstruct.valdiff = valdiffunique;
outputstruct.valdiff_abs = valdiffunique_abs;
outputstruct.item1chosen_valdiff_all = item1chosen_valdiff_all;
outputstruct.item1chosen_valdiff_lastL = item1chosen_valdiff_lastL;
outputstruct.item1chosen_valdiff_lastR = item1chosen_valdiff_lastR;
outputstruct.item1chosen_timeadv_norm = item1chosen_timeadv_norm;
outputstruct.item1chosen_timeadv_beta = item1chosen_timeadv_beta;
outputstruct.firstchosen_firstfixdur = firstchosen_firstfixdur;
outputstruct.firstchosen_firstfixdur_norm = firstchosen_firstfixdur_norm;
outputstruct.fixdiff_valdiff = fixdiff_valdiff;

outputstruct.rt = rt_mean;
outputstruct.rt_valdiff_abs = rt_valdiff_abs;
outputstruct.switchcount_valdiff_abs = switchcount_valdiff_abs;
outputstruct.switchcountRT_valdiff_abs = switchcountRT_valdiff_abs;

outputstruct.binmeans_item1timeadv = means_item1timeadv;
outputstruct.binmeans_firstfixdur = means_firstfixdur;

outputstruct.mid_fixitemval = mid_fixitemval;
outputstruct.mid_fixunfixdiff = mid_fixunfixdiff;
outputstruct.mid_absvaldiff = mid_absvaldiff;



end

