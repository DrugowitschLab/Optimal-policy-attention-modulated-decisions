function dstruct = getEmptyDstruct_realData(dstruct_real,iter)
%Get empty dstruct with same subjects, trials and item values as inputted dstruct
% iter: repeat the same trials x times to increase trial number

dstruct = struct;

for s = 1:length(dstruct_real)
    numtrials = length(dstruct_real(s).trialnum) * iter;
    for ff = fieldnames(dstruct_real)'
        if iscell(dstruct_real(1).(char(ff)))
            dstruct(s).(char(ff)) = cell(numtrials,size(dstruct_real(1).(char(ff)),2));
        elseif isnumeric(dstruct_real(1).(char(ff)))
            dstruct(s).(char(ff)) = nan(numtrials,size(dstruct_real(1).(char(ff)),2));
        end
    end
    dstruct(s).trialnum = [1:numtrials]';
    dstruct(s).itemval = repmat(dstruct_real(s).itemval,iter,1);
end

end

