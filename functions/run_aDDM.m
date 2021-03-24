function [choice,rt,fixitem,fixdur,tItem,yseq] = run_aDDM(z,g,dt,k,sig2,decbound,yseq_toMaxTime)
%Compute the relative decision value as shown in Krajbich et al., 2010
% Inputs:
% yseq = sequence of fixations (1,2) to the dec time limit. This will be
% trimmed in the output based on the response time

% Paper nomenclature: 
% Viewing right: V_t = V_{t-1} + d*(r_left - θ*r_right) + e_t
% Viewing left: V_t = V_{t-1} + d*(-θ*r_left + r_right) + e_t
% e_t = N(0,σ^2)
% Best fitting model: d = 0.0002ms−1 (d = 0.2s-1), θ = 0.3 and σ = 0.02

% Our notation:
% d = k*dt
% σ^2 = sig2*dt

rdv = 0;
decmade = false;
for ni = 1:length(yseq_toMaxTime)
    y = yseq_toMaxTime(ni);
    rdv = rdv + (k*dt)*( z(1)*(g^(y-1)) - z(2)*(g^(2-y)) ) + randn*sqrt(sig2*dt);
    if abs(rdv) >= decbound
        decmade = true;
        break;
    end
end
yseq = yseq_toMaxTime(1:ni);

% get behavior
if decmade, choice = double(rdv<0) + 1;
else, choice = nan;
end
rt = ni*dt;
fixitem = yseq([1,diff(yseq)]~=0);
fixdur = diff([find([1,diff(yseq)]~=0),ni+1])*dt;
tItem = [sum(fixdur(fixitem==1)),sum(fixdur(fixitem==2))];

if any(isnan(fixitem)), keyboard; end

end

