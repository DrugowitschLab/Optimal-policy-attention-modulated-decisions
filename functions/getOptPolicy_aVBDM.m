function optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir)
%Computes the value function using Bellman's equation for attention-based
%value-based decision making
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
    fprintf('\nComputing optimal policy: %s\n',saveStrAppend);
    N = length(params.ts);
    
    % Compute the value function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get value for deciding (Vd)
    R1 = nan(length(params.zs),length(params.zs));  % Vd for y=1
    R2 = nan(length(params.zs),length(params.zs));  % Vd for y=2
    for z1i = 1:length(params.zs)
        for z2i = 1:length(params.zs)
            R1(z1i,z2i) = round(( params.zs(z1i)-params.zs(z2i) )/2,10);
            R2(z1i,z2i) = round(( params.zs(z2i)-params.zs(z1i) )/2,10);
        end
    end
    bothparams.Rs = {R1,R2};
    
    % Show the relationship between the item values and the Delta term
    makeplots_Vd = 0;
    if makeplots_Vd==1
        % Show how the Vd is determined for both y=1 (1st row) and y=2 (2nd row)
        colors = {'autumn','winter'};
        figure; spIdx = [1,5];
        for y = 1:2
            subplot(2,4,spIdx(y));
            surf(params.zs,params.zs,bothparams.Rs{y}); shading flat; axis xy; hold on; pbaspect([1,1,1]);
            set(gca,'zlim',[min(params.Rs),max(params.Rs)]);
            colormap(colors{y}); freezeColors;
            title('r_1-r_2'); xlabel('z2'); ylabel('z1'); zlabel('Delta');
            subplot(2,4,spIdx(y)+1);
            surf(params.zs,params.zs,bothparams.Rs{3-y}); shading flat; axis xy; hold on; pbaspect([1,1,1]);
            set(gca,'zlim',[min(params.Rs),max(params.Rs)]);
            colormap(colors{3-y}); freezeColors;
            title('r_2-r_1'); xlabel('z2'); ylabel('z1'); zlabel('Delta');
            subplot(2,4,spIdx(y)+2);
            surf(params.zs,params.zs,bothparams.Rs{y}); shading flat; axis xy; hold on;
            set(gca,'zlim',[min(params.Rs),max(params.Rs)]);
            colormap(colors{y}); freezeColors;
            surf(params.zs,params.zs,bothparams.Rs{3-y}); shading flat;
            colormap(colors{3-y}); freezeColors;
            title('Both'); xlabel('z2'); ylabel('z1');  zlabel('Delta'); pbaspect([1,1,1]);
            subplot(2,4,spIdx(y)+3);
            surf(params.zs,params.zs,max(bothparams.Rs{y},bothparams.Rs{3-y})); shading flat; axis xy; hold on;
            set(gca,'zlim',[min(params.Rs),max(params.Rs)]);
            title('Vd'); xlabel('z2'); ylabel('z1');  zlabel('Delta'); pbaspect([1,1,1]);
            colormap(linspecer);
        end
    end
    
    % Get the decision boundary for deciding immediately (Vd) for both y conditions
    % Take the diagonal cut of the planes to reduce to lower dimension, along (r1+r2)/2
    Vd = [];
    for y = 1:2
        Vd_lowdim = [];
        Vd_y = max(bothparams.Rs{y},bothparams.Rs{3-y});
        for di = size(Vd_y,1):-1:-size(Vd_y,1)
            Vd_lowdim = cat(1,Vd_lowdim,unique(diag(Vd_y,di)));
        end
        Vd = cat(1,Vd,Vd_lowdim');
    end
    
    val = struct;
    val.Vd_1 = Vd(1,:); val.Vd_2 = Vd(2,:);
    val.Ve_1 = nan(N,N,length(params.Rs)); val.Ve_2 = nan(N,N,length(params.Rs));  % Value for  deciding later (accumulating more evidence)
    val.V1 = nan(N,N,length(params.Rs)); val.V2 = nan(N,N,length(params.Rs));  % Value for deciding
    val.p1 = nan(N,N,length(params.Rs)); val.p2 = nan(N,N,length(params.Rs));  % Optimal policy (1,2:choose options [1,2], 3: accumulate, 4: switch)
    
    % Above certain tMax = t1+t2, Value = Vd
    n1n2 = [nchoosek(1:N,2);repmat([1:N]',1,2)];  % x and y indices
    n1n2 = n1n2(sum(n1n2,2) >= N,:);
    n1n2 = unique([n1n2;fliplr(n1n2)],'rows');
    for i = 1:size(n1n2,1)
        for y = 1:2
            val.(sprintf('V%d',y))(n1n2(i,1),n1n2(i,2),:) = val.(sprintf('Vd_%d',y));
            [~,Vd_min_i] = min(val.(sprintf('Vd_%d',y)));
            val.(sprintf('p%d',y))(n1n2(i,1),n1n2(i,2),1:Vd_min_i-1) = 2;
            val.(sprintf('p%d',y))(n1n2(i,1),n1n2(i,2),Vd_min_i+1:end) = 1;
        end
    end
    
    % Compute value function
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
    every10Prc = floor(prctile(1:size(n1n2,1),10:10:100));
    for i = 1:size(n1n2,1)
        if any(i==every10Prc), fprintf('%.0f%% ',find(i==every10Prc)*10); end
        n1 = n1n2(i,1); n2 = n1n2(i,2);
        t1 = n1*params.dt; t2 = n2*params.dt;
        
        % Transition matrix: T = p(R(t+dt)|R(t),t1,t2,y)
        % y-axis: R(t+dt), x-axis: R(t)
        tMat = get_tMat(params.Rs,t1,t2,params);
        
        % Compute the value of accumulating evidence
        Ve_1 = squeeze(val.V1(n1+1,n2,:))'*tMat(:,:,1) - params.c*params.dt;
        Ve_2 = squeeze(val.V2(n1,n2+1,:))'*tMat(:,:,2) - params.c*params.dt;
        % Compute the values for V1 and V2, given the value for accumulating
        % evidence and cost for switching
        [val.V1(n1,n2,:),val.V2(n1,n2,:),val.p1(n1,n2,:),val.p2(n1,n2,:)] = computeValue(val.Vd_1,val.Vd_2,Ve_1,Ve_2,params.Rs,params.Cs);
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
    n1=20; n2 = 20;
    if ~isequal_aij(squeeze(val.V1(n1,n2,:)),squeeze(val.V2(n2,n1,:)),1)
        warning('V1 and V2 are NOT symmetric');
    end
    
    % only output the optimal policy
    optpolicy = struct;
    optpolicy.p1 = val.p1;
    optpolicy.p2 = val.p2;
    
    save(fullfile(optpolicy_dir,datafilename),'optpolicy','params','-v7.3');  % Save optimal policy
end

end

