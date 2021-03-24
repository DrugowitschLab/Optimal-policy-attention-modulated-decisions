function [choice,rt,tItem,fixseq,fixdur,fixitem,reward,delta_opt] = runTrial_aVBDM(z,optpolicy,params)
%Run a single trial in the attention value-based decision making paradigm
% input:
% val: structure with policy information
% params: basic parameters of task
% params_fb: optional input to simulate with flat boundary
% output:
% choice: option 1 or 2
% rt: response time
% valGain: total value gained from trial
% fixseq: sequence of attended items

N = length(params.ts);

% If z is not provided, generate a random z using the parameters
if isnan(z)
    output = generateTrial(params,params.zbar);
    z = output.z;
end

choice = double(rand>0.5)+1;  % default is random choice (in case model ends in policy space where delta = 0)

dx = nan(N,2);   % Evidence accumulation
mean_post = nan(N,2);   % Posterior mean
delta_opt = nan(N,1);  % Delta ((r1 - r1) / 2)
tItem = [0,0];   % Time points each item was attended to
fixseq = nan(1,N);   % Sequence of attended items
y = (rand>0.5)+1;   % Currently attended item
n = 1; decmade = false;
while n < N && ~decmade
    tItem(y) = tItem(y)+params.dt;  % Time advances
    fixseq(n) = y;
    [dx,mean_post(n,:)] = getPosterior_aVBDM(n,dx,fixseq,tItem,z,params);
    % subtract cost of choosing nonfixated item
    currenttrialpost = mean_post(n,:);
    delta_opt(n) = (currenttrialpost(1)-currenttrialpost(2))/2;
    
    % Policy: 1,2 (choose option 1,2), 3 (accumulate), 4 (switch)
    t1_i = isequal_aij(params.ts,tItem(1),2); t2_i = isequal_aij(params.ts,tItem(2),2);  % Deal with precision error
    [~,R_i] = min(abs(params.Rs-delta_opt(n)));
    trialOptPolicy = optpolicy.(sprintf('p%d',y))(t1_i,t2_i,R_i);
    if trialOptPolicy==1 || trialOptPolicy==2
        choice = trialOptPolicy;
        decmade = true;
    else
        if trialOptPolicy==4  % Switch attention
            y = 3-y;
        end
        n = n+1;
    end
end
rt = n*params.dt;
fixseq = fixseq(~isnan(fixseq));

% organize fixation properties
switchPoints = find(diff(fixseq(~isnan(fixseq)))~=0)+1;
switchPointsMod = [1,switchPoints,sum(~isnan(fixseq))+1];
% Get the fixated items and the fixation duration
fixdur = nan(1,length(switchPoints)+1); fixitem = nan(1,length(switchPoints)+1);
for sp_i = 1:length(switchPointsMod)-1
    fixdur(sp_i) = (switchPointsMod(sp_i+1)-switchPointsMod(sp_i))*params.dt;
    fixitem(sp_i) = fixseq(switchPointsMod(sp_i));
end

reward = z(choice) - length(switchPoints)*params.Cs - rt*params.c;

end

