function optpolicy = getOptPolicy_aPDM_R(params,optpolicy_dir)
%Computes the value function using Bellman's equation for attention-based
%perceptual decision making
% Inputs:
% params: task parameters
% optpolicy_dir: directory where optimal policy is saved, in case this was
% already computed. If not, the function will save the computed policy
% space into this directory
% Outputs:
% 1. val: state-space with values for each state
% 2. optimal boundaries for behavior (choice, accumulate, switch)


% Check if policy space is already computed and saved. If so, just load.
saveStrAppend = getOptPolicySaveStrAppend(params);
datafilename = sprintf('optPolicyData%s.mat',saveStrAppend);
alldatafiles = getDirFileNames(optpolicy_dir);

% First, check if data is already stored
if ismember(datafilename,alldatafiles)
    load(fullfile(optpolicy_dir,datafilename));
else
    N = length(params.ts);
    s = params.sig2/params.sig2_z;
    
    % Initialize structure to store optimal policy data
    nanMat = nan(N,N,length(params.Rs));
    val = struct;
    val.Vd = nanMat;
    val.Ve_1 = nanMat; val.Ve_2 = nanMat;
    val.p1 = nanMat; val.p2 = nanMat;
    val.V1 = nanMat; val.V2 = nanMat;
    
    % Compute the value function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute value for deciding for each R,t1,t2 combination
    for n1 = 1:N
        for n2 = 1:N
            sig2_1 = params.sig2/(s+(n1*params.dt)+params.aGamma*(n2*params.dt));
            sig2_2 = params.sig2/(s+params.aGamma*(n1*params.dt)+(n2*params.dt));
            Vd_1 = normcdf((2*params.Rs)./sqrt(sig2_1+sig2_2));
            val.Vd(n1,n2,:) = max(Vd_1,1-Vd_1);
            
            % Above certain tMax = t1+t2, Value = Vd
            if n1+n2 >= N
                for y = 1:2
                    val.(sprintf('V%d',y))(n1,n2,:) = val.Vd(n1,n2,:);
                    [~,Vd_min_i] = min(val.Vd(n1,n2,:));
                    val.(sprintf('p%d',y))(n1,n2,1:Vd_min_i-1) = 2;
                    val.(sprintf('p%d',y))(n1,n2,Vd_min_i+1:end) = 1;
                end
            end
        end
    end
    
    % Compute value function diagonally, so that V(t1+t2 = tmax) = Vd for both V1,V2
    % Find indices for each diagonal
    idxMat = reshape(1:N^2,N,N);
    diagIdx = nan(N,N,N-12); diagCounter = 1;
    for diagNumber = 2:N-1
        diagIdx(:,:,diagCounter) = fliplr(ismember(idxMat,diag(idxMat,-diagNumber))');
        diagCounter = diagCounter+1;
    end
    n1n2 = [];
    for d = 1:size(diagIdx,3)
        [row,col]=find(diagIdx(:,:,d)==1);
        n1n2 = cat(1,n1n2,[row,col]);
    end
    
    % Use backward induction to compute the value function for all t1, t2,
    % and Delta
    % Starting from V(t1+t2 = tmax) = Vd for both V1,V2, go down diagonally
    % and compute the value function
    fprintf('\nComputing value function diagonally: ');
    every10Prc = floor(prctile(1:size(n1n2,1),10:10:100));
    for i = 1:size(n1n2,1)
        if any(i==every10Prc), fprintf('%.0f%% ',find(i==every10Prc)*10); end
        n1 = n1n2(i,1); n2 = n1n2(i,2);
        t1 = n1*params.dt; t2 = n2*params.dt;
        
        % Transition matrix: T = p(R(t+dt)|R(t),t1,t2,y)
        % y-axis: R(t+dt), x-axis: R(t)
        tMat = get_tMat(params.Rs,t1,t2,params);
        
        % value for deciding
        Vd = squeeze(val.Vd(n1,n2,:))';
        % Compute the value of accumulating evidence
        Ve_1 = squeeze(val.V1(n1+1,n2,:))'*tMat(:,:,1) - params.c*params.dt;
        Ve_2 = squeeze(val.V2(n1,n2+1,:))'*tMat(:,:,2) - params.c*params.dt;
        % Compute the values for V1 and V2, given the value for accumulating
        % evidence and cost for switching
        [val.V1(n1,n2,:),val.V2(n1,n2,:),val.p1(n1,n2,:),val.p2(n1,n2,:)] = computeValue(Vd,Vd,Ve_1,Ve_2,params.Rs,params.Cs);
        val.Ve_1(n1,n2,:) = Ve_1;
        val.Ve_2(n1,n2,:) = Ve_2;
    end
    fprintf('\n');
    clear Ve_1 Ve_2;
    
    % Validate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % V1-V2 should be bounded by Cs
    temp = val.V1-val.V2;
    [lowBound,highBound] = bounds(temp(:));
    if abs(lowBound)-params.Cs > 1.0e-8 || abs(highBound)-params.Cs > 1.0e-8
        warning('V1-V2 is not bounded by Cs');
    end
    clear temp;
    
    % Check for symmetry
    % V(R,t1,t2,1) = V(R,t2,t1,2)
    n1=75; n2 = 50;
    if ~isequal_aij(squeeze(val.V1(n1,n2,:)),squeeze(val.V2(n2,n1,:)),1)
        warning('V1 and V2 are NOT symmetric');
    end
    
    % Just save the optimal policy
    optpolicy = struct;
    optpolicy.p1 = val.p1;
    optpolicy.p2 = val.p2;
    
    % Save optimal policy that can be loaded next time
%     save(fullfile(optpolicy_dir,datafilename),'optpolicy','params','-v7.3');  % Save optimal policy
end
