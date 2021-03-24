function dstruct = getEmptyDstruct(dstruct_real,numsubs,z_all)
%Get empty dstruct with specified number of subjects and trials

dstruct = struct;
dstruct_sub = struct;
numtrials = size(z_all,1);
for ff = fieldnames(dstruct_real)'
    if iscell(dstruct_real(1).(char(ff)))
        dstruct_sub.(char(ff)) = cell(numtrials,size(dstruct_real(1).(char(ff)),2));
    elseif isnumeric(dstruct_real(1).(char(ff)))
        dstruct_sub.(char(ff)) = nan(numtrials,size(dstruct_real(1).(char(ff)),2));
    end
end
dstruct_sub.itemval = z_all;

% New field: delta for showing particle trajectory
dstruct_sub.delta_opt = cell(numtrials,1);

for s = 1:numsubs
    if s==1, dstruct = dstruct_sub;
    else, dstruct(s) = dstruct_sub;
    end
end



end

