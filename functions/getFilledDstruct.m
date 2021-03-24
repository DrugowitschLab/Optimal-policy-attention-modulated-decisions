function dstruct = getFilledDstruct(dstruct_real,numsubs,z_all)
%Get filled dstruct with specified number of subjects. 
% uses input of z pair values, and samples fixation behavior randomly from
% data

% Compile data into large array to easily sample randomly for each
% simulated trial
numtrials = size(z_all,1);
fixdur_all = {};
fixitem_all = {};
tItem_all = [];
rt_all = [];
for s = 1:length(dstruct_real)
    fixdur_all = cat(1,fixdur_all,dstruct_real(s).fixdur);
    fixitem_all = cat(1,fixitem_all,dstruct_real(s).fixitem);
    tItem_all = cat(1,tItem_all,dstruct_real(s).tItem);
    rt_all = cat(1,rt_all,dstruct_real(s).rt);
end

% Initialize empry struct
dstruct = struct;
dstruct_sub = struct;
for ff = fieldnames(dstruct_real)'
    if iscell(dstruct_real(1).(char(ff)))
        dstruct_sub.(char(ff)) = cell(numtrials,size(dstruct_real(1).(char(ff)),2));
    elseif isnumeric(dstruct_real(1).(char(ff)))
        dstruct_sub.(char(ff)) = nan(numtrials,size(dstruct_real(1).(char(ff)),2));
    end
end
for s = 1:numsubs
    if s==1, dstruct = dstruct_sub;
    else, dstruct(s) = dstruct_sub;
    end
end

% Fill struct
for s = 1:numsubs
    dstruct(s).trialnum = [1:numtrials]';
    dstruct(s).itemval = z_all;
    dstruct(s).choice = nan(numtrials,1);
    for ti = 1:size(z_all,1)
        i_rand = randperm(size(rt_all,1),1);
        dstruct(s).trialnum(ti) = ti;
        dstruct(s).fixdur{ti} = fixdur_all{i_rand};
        dstruct(s).fixitem{ti} = fixitem_all{i_rand};
        dstruct(s).rt(ti) = rt_all(i_rand);
        dstruct(s).tItem(ti,:) = tItem_all(i_rand,:);
    end
end

end

