function perf = getPerfFromDstruct(dstruct,params)
%From data structure dstruct, return the per-subject performance measures
perf = struct;
perf.acc = nan(length(dstruct),1);
perf.rt = nan(length(dstruct),1);
perf.reward = nan(length(dstruct),1);
for s = 1:length(dstruct)
    % Ignore trials where z1 = z2?
%     i_trials_use = dstruct(s).itemval(:,2)-dstruct(s).itemval(:,1) ~= 0;  % ignore z1=z2 trials
    i_trials_use = true(size(dstruct(s).itemval,1),1);
    % Accuracy
    choice_corr = double(dstruct(s).itemval(i_trials_use,2)-dstruct(s).itemval(i_trials_use,1) > 0) + 1;  % correct choice
    choice_data = dstruct(s).choice(i_trials_use);
    % Mean single-trial reward
    numswitches = cellfun('length',dstruct(s).fixdur(i_trials_use))-1;
    idx_choice = [dstruct(s).choice==1,dstruct(s).choice==2];
    itemval_chosen = sum(dstruct(s).itemval.*idx_choice,2);
    rt = dstruct(s).rt(i_trials_use);
    reward_all = itemval_chosen(i_trials_use,:) - (rt.*params.c + numswitches.*params.Cs);
    
    % compile data
    if length(dstruct) > 1
        perf.acc(s) = mean(choice_corr==choice_data);  % Accuracy
        perf.rt(s) = mean(rt);  % RT
        perf.reward(s) = mean(reward_all);  % Mean single-trial reward
    else
        % if it's a single subject, just return all the trial information
        perf.acc = choice_corr==choice_data;  % Accuracy
        perf.rt = rt;  % RT
        perf.reward = reward_all;  % Mean single-trial reward
    end
end

end

