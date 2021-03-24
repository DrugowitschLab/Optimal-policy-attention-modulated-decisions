function [V1,V2,p1,p2] = computeValue(Vd_1,Vd_2,Ve_1,Ve_2,Rs,Cs)
%Returns the value and policy for options 1 and 2 given the Vd, Ve, and cost for
%switching (C_s)
%Cd: cost for choosing item not attended to
% Policy: 1,2 (choose options 1,2), 3 (accumulate), 4 (switch)

% None of the values should be nans
if any(isnan(Vd_1)) || any(isnan(Vd_2)) || any(isnan(Ve_1)) || any(isnan(Ve_2))
    error('There are NaNs');
end

% Initially:
% 1: choice
% 2: accumulate
% 3: switch
% Compute V1, assume no switch
[V1,p1]=max([Vd_1; Ve_1],[],1);

% Compute V2
[V2,p2] = max([Vd_2; Ve_2; V1-Cs],[],1);

% Check if V2(t) was determined by V1(t)-C_s
recompute_i = abs(V2-(V1-Cs)) > 1.0e-8;
if any(recompute_i)
    [V1(recompute_i),p1(recompute_i)] = max([Vd_1(recompute_i); Ve_1(recompute_i); V2(recompute_i)-Cs],[],1);
end


% Convert policy so that:
% 1,2 --> choose options 1,2
% 3 --> accumulate
% 4 --> switch
p1 = p1+1; p2 = p2+1;

% Rs > 0: R(1) > R(2), vice versa
Rs_choose1 = Rs'>0;
p1(p1==2 & Rs_choose1) = 1;
p2(p2==2 & Rs_choose1) = 1;

end

