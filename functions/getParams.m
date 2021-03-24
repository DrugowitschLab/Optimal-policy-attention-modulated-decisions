function params = getParams(whichTask,paramInput)
% Parameters for value-based decision making
% zbar = prior mean
% switchRate = Chance of switching attention at each time point
% Presence of switchRate input simulates attentional switches

params = struct;
params.zsRes = 0.05;
params.zbar = 0;
params.dt = 0.05;
params.t_s = params.dt;
params.T = 6;
params.ts = 0:params.dt:params.T;
if nargin==2  % input parameter values
    params = paramInput;
    
elseif whichTask==1
    params.zs = -10:params.zsRes:10;
    params.Rs = unique(round(diff(permn(params.zs,2),1,2)/2,10));

    % PAPER PARAMS
    % sortmethod = 3; mineffect = [5,8,1];
    params.sig2 = 27;
    params.sig2_z = 18;
    params.c = 0.23;
    params.Cs = 0.018;
    params.c_s = params.Cs - params.t_s*params.c;
    params.aGamma = 0.004;

elseif whichTask==2
    params.zs = -10:params.zsRes:10;
    params.Rs = unique(round(diff(permn(params.zs,2),1,2)/2,10));
    
    % sortmethod = 3; mineffect = [10,10,8];
    params.sig2 = 21;
    params.sig2_z = 54;
    params.c = 0.06;
    params.Cs = 0.013;
    params.c_s = params.Cs - params.t_s*params.c;
    params.aGamma = 0.009;
end

% Covariance matrices - Assume two variables are uncorrelated
params.cMatZ = diag([params.sig2_z,params.sig2_z]);
params.cMat = diag([params.sig2,params.sig2]);



end
