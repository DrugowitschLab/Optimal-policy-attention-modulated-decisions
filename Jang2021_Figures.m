%% Creat main figures for Jang et al., 2021
% Optimal policy for attention-modulated decisions explains human fixation behavior
% Run this block first, then any other block below to make the figures
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET CODE DIRECTORY HERE TO RUN ALL SCRIPTS
rootdir = '/Users/anthony_jang/Dropbox (MIT)/Research/Ongoing/att_paper/202009-elife_submission/Code';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(rootdir));
data_dir = fullfile(rootdir,'data');
optpolicy_dir = fullfile(rootdir,'data/optimal_policy_save');
set(0,'defaultAxesFontSize',20);


%% Figure 1C
%%%%%%%%%%%%
whichTask = 1;  % 1: Value-based task, 2: Perceptual task
% Load params
params = getParams(whichTask);
% Parameters were chosen empirically to best show the evidence accumulation
% process.
z = [13,10];  % Item values
sig2 = 5;  % Evidence accumulation noise
sig2_z = 10 ;  % Prior variance
zbar = 0;  % Mean of the prior
aGamma = 0.1;  % Inattention noise
dt = 0.01;  % Time step
t = 3;  % Length of time to show
ts = dt:dt:t;
N = t/dt;
colors_item = linspecer(4);
colors_item = colors_item([1,2],:);  % blue and green items

fixseq = [2*ones(1,round((4*N)/12)),ones(1,round((3*N)/12)),2*ones(1,round((3*N)/12))];
fixseq = [fixseq,ones(1,N-length(fixseq))];
switchpts = ts([0,diff(fixseq)~=0]==1);

% Plot evidence & posterior distribution across time
dx = nan(N,2);   % Evidence accumulation
post_mean = nan(N,2);   % Posterior mean
post_sig2 = nan(N,2);   % Posterior variance
tItem = [0,0];   % Time points each item was attended to
for n = 1:N
    y = fixseq(n);  % Currently attended item, [1,2]
    tItem(y) = tItem(y)+dt;  % Time advances
    % Sequence of attntion weights per time point
    aGammaSeq_1 = aGamma.*(fixseq(1:n)==2) + double(fixseq(1:n)==1);
    aGammaSeq_2 = aGamma.*(fixseq(1:n)==1) + double(fixseq(1:n)==2);
    % Evidence accumulation
    dx1 = randn*sqrt( (sig2/(aGamma^(y-1)))*dt ) + z(1)*dt;
    dx2 = randn*sqrt( (sig2/(aGamma^(2-y)))*dt ) + z(2)*dt;
    dx(n,:) = [dx1,dx2];
    % Mean of the posterior distribution
    X1 = sum(dx(1:n,1).*aGammaSeq_1');
    post_mean(n,1) = (zbar*sig2 + X1*sig2_z) / (sig2 + (tItem(1)+aGamma*tItem(2))*sig2_z);
    post_sig2(n,1) = sig2/(sig2/sig2_z + tItem(1)+aGamma*tItem(2));
    X2 = sum(dx(1:n,2).*aGammaSeq_2');
    post_mean(n,2) = (zbar*sig2 + X2*sig2_z) / (sig2 + (aGamma*tItem(1)+tItem(2))*sig2_z);
    post_sig2(n,2) = sig2/(sig2/sig2_z + aGamma*tItem(1)+tItem(2));
end
% Plot
fh = figure('units','normalized','outerposition',[1,1,0.5,1]);
% Show momentary evidence across time
plots = [];
subplot(3,1,1); hold on;
for i = 1:2
    if i==1, plots(i) = plot(ts,dx(:,i),'o','Color',colors_item(i,:),'markersize',6);
    else, plots(i) = plot(ts,dx(:,i),'.','Color',colors_item(i,:),'markersize',12);
    end
    plot(get(gca,'xlim'),[z(i)*dt,z(i)*dt],'--','Color',colors_item(i,:));
end
set(gca,'ylim',[-6,6]);
for sw = 1:length(switchpts), plot([switchpts(sw),switchpts(sw)],get(gca,'ylim'),'k--'); end
ylabel('Momentary evidence');
legend(plots,{'Item 1','Item 2'});
% Show posterior across time
plots = [];
subplot(3,1,2); hold on;
for i = 1:2
    plots(i) = shadedErrorBars(ts,post_mean(:,i)',sqrt(post_sig2(:,i))',colors_item(i,:),0.5,'-');
    plot(get(gca,'xlim'),[z(i),z(i)],'--','Color',colors_item(i,:));
end
% set(gca,'ylim',[-10,20]);
for sw = 1:length(switchpts), plot([switchpts(sw),switchpts(sw)],get(gca,'ylim'),'k--'); end
xlabel('Time (s)'); ylabel('Posterior');
legend(plots,{'Item 1','Item 2'});
% Show posterior variance across time
plots = [];
subplot(3,1,3); hold on;
for i = 1:2
    plots(i) = plot(ts,post_sig2(:,i),'Color',colors_item(i,:));
end
for sw = 1:length(switchpts), plot([switchpts(sw),switchpts(sw)],get(gca,'ylim'),'k--'); end
legend(plots,{'Item 1','Item 2'});
xlabel('Time (s)'); ylabel('Posterior sd');



%% Figure 2
%%%%%%%%%%%
whichTask = 1;  % 1: Value-based, 2: Perceptual
params = getParams(whichTask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the optimal policy in 3D space and 2D slices

% Get optimal policy
optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir);

% Plot settings
zlimit = [-3,3];
zval = params.Rs;
zStr = '\Delta';
colors_policy = linspecer(4); colors_y = [1,0,0;0,0,1];
% Get boundaries for optimal policy
optbound = getOptpolicyBoundary(optpolicy,params);
xlimit = [min(params.ts),max(params.ts)]; ylimit = [min(params.ts),max(params.ts)];
deltastr = '\Delta';

y = 1;  % Show the policy space for y = 1 (attending to item 1)
R_slices = [-1,0,1];  % Cross-sections to show

% Plot slices
spIdx = 1;
figure;
for R = R_slices
    R_i = zval==R;
    thisslice = optpolicy.(sprintf('p%d',y))(:,:,R_i);
    if R==0
        % 1. Red & blue stripes version
        stripeNum = 6;
        stripeIdx = discretize(1:size(thisslice,1),[-inf,quantile(1:size(thisslice,1),stripeNum-1),inf]);
        stripeIdxMat = zeros(size(thisslice));
        stripeIdxMat_odd = stripeIdxMat;
        stripeIdxMat_odd(:,mod(stripeIdx,2)==1) = 1;
        thisslice(stripeIdxMat_odd & isnan(thisslice)) = 1;
        thisslice(~stripeIdxMat_odd & isnan(thisslice)) = 2;
    end
    subplot(1,3,spIdx); imagesc(params.ts,params.ts,thisslice); hold on; axis xy;
    policyColorSetting;
    title(sprintf('Delta(t) = %d',R)); xlabel('t2'); ylabel('t1');
    pbaspect([1,1,1]);
    spIdx = spIdx+1;
end

% Plot 3D optimal policy boundaries
y = 2;
axison = 1;
sliceon = 1;  % Show gray squares that "slice" the policy space
figure('units','normalized','outerposition',[0 0 0.4 0.7]);
fNames = {'dec','accum','switch'};
for fi = 2:length(fNames)
    if strcmp(fNames{fi},'accum'), policycolor = colors_policy(4,:);
    elseif strcmp(fNames{fi},'switch'), policycolor = colors_policy(3,:);
    end
    if ~isempty(optbound.(fNames{fi}){y}) % a policy might not exist (e.g. no switch boundaries)
        R_i = ~isnan(optbound.(fNames{fi}){y}(:,3));
        % Remove duplicate points
        xx = params.ts(optbound.(fNames{fi}){y}(R_i,1))';
        yy = params.ts(optbound.(fNames{fi}){y}(R_i,2))';
        zz = optbound.(fNames{fi}){y}(R_i,3);
        hold on;
        
        % Run  program
        trisurface=MyCrustOpen([xx,yy,zz]);
        
        % plot the output triangulation
        trisurf(trisurface,xx,yy,zz,'facecolor',policycolor,'edgecolor','none','facealpha',1);
        grid on;
        if axison==0, axis off; end
        light_handle = light;
        lightangle(light_handle,-30,0);
        material dull;
        set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit);
        ylabel('t1[s]'); xlabel('t2[s]'); zlabel(sprintf('%s(t)',deltastr))
        pbaspect([range(xlimit),range(ylimit),range(zlimit)]);
    end
end
set(gca,'xtick',0:2:6,'ytick',0:2:6,'ztick',-2:2);
setFigFontSize(26,fh);
lightangle(light_handle,80,0);

if sliceon==1
    [X,Y] = meshgrid(params.ts,params.ts);
    for zz = 1:length(R_slices)
        s = surface(X,Y,ones(size(X)).*R_slices(zz));
        s.FaceColor = 'k';
        s.FaceAlpha = 0.3;
        s.EdgeColor = 'none';
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show particle trajectory for an example trial
% Load optimal policy information
optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir);

% Simulate behavior. Iterate until a trial with 3+ fixations
foundgoodtrial = false;
while ~foundgoodtrial
    z = [1,1];
    numtrials = 100;
    particle_xyz = [0,0,0];
    N = length(params.ts);
    dx = nan(N,2);   % Evidence accumulation
    mean_post = nan(N,2);   % Posterior mean
    tItem = [0,0];   % Time points each item was attended to
    tItem_all = [0,0];
    fixseq = nan(1,N);   % Sequence of attended items
    y = 1;
    n = 1; decmade = false;
    fixnum = nan(1,N); fn = 1;
    while n <= N && ~decmade
        tItem(y) = tItem(y)+params.dt;  % Time advances
        tItem_all = cat(1,tItem_all,tItem);
        fixseq(n) = y;
        fixnum(n) = fn;
        [dx,mean_post(n,:)] = getPosterior_aVBDM(n,dx,fixseq,tItem,z,params);
        % subtract cost of choosing nonfixated item
        currenttrialpost = mean_post(n,:);
        R_mean = (currenttrialpost(1)-currenttrialpost(2))/2;
        
        % Policy: 1,2 (choose option 1,2), 3 (accumulate), 4 (switch)
        t1_i = isequal_aij(params.ts,tItem(1),2); t2_i = isequal_aij(params.ts,tItem(2),2);  % Deal with precision error
        [~,R_i] = min(abs(params.Rs-R_mean));
        trialpolicy = optpolicy.(sprintf('p%d',y))(t1_i,t2_i,R_i);
        if trialpolicy==1 || trialpolicy==2, decmade = true;
        else
            if trialpolicy==4  % Switch attention
                y = 3-y;
                fn = fn+1;
            end
            n = n+1;
        end
        particle_xyz = cat(1,particle_xyz,[tItem(2),tItem(1),R_mean]);
    end
    
    if max(fixnum)==3, foundgoodtrial = true; end
end

% Plot reconstructed surface
colors_policy = linspecer(4); colors_y = [1,0,0;0,0,1];
xlimit = [-1,max(params.ts)]; ylimit = [-1,max(params.ts)]; zlimit = [-3,3];
optbound = getOptpolicyBoundary(optpolicy,params);
y = 1;
axison = 1;
deltastr = '\Delta';
fixnumunique = unique(fixnum(~isnan(fixnum)));
if length(fixnumunique)>3, fixnumunique = fixnumunique(1:3); end
for fn = fixnumunique
    fh1 = figure; hold on;
    
    fNames = {'dec','accum','switch'};
    for fi = 2:length(fNames)
        if strcmp(fNames{fi},'accum'), policycolor = colors_policy(4,:);
        elseif strcmp(fNames{fi},'switch'), policycolor = colors_policy(3,:);
        end
        if ~isempty(optbound.(fNames{fi}){y}) % a policy might not exist (e.g. no switch boundaries)
            R_i = ~isnan(optbound.(fNames{fi}){y}(:,3));
            % Remove duplicate points
            xx = params.ts(optbound.(fNames{fi}){y}(R_i,1))';
            yy = params.ts(optbound.(fNames{fi}){y}(R_i,2))';
            zz = optbound.(fNames{fi}){y}(R_i,3);
            
            % Run program to reconstruct 3D surface
            trisurface=MyCrustOpen([xx,yy,zz]);
            
            % Plot the output triangulation
            trisurf(trisurface,xx,yy,zz,'facecolor',policycolor,'edgecolor','none','facealpha',1);
            grid on;
            if axison==0, axis off; end
            light_handle = light;
            lightangle(light_handle,-30,0);
            material dull;
            set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit);
            ylabel('t1[s]'); xlabel('t2[s]'); zlabel(sprintf('%s(t)',deltastr));
            pbaspect([1,1,1]);
        end
    end
    lightangle(light_handle,80,0);
    
    i_fixnum = find(fixnum==fn);
    % End of this fixation overlabs with beginning of next fixation.
    % Take the first timestep of the next fixation
    if i_fixnum(end) < length(fixnum), i_fixnum = [i_fixnum,i_fixnum(end)+1]; end
    
    if fn==1, t_fixstart = 0;
    else, t_fixstart = tItem_all(i_fixnum(1),3-y);
    end
    if y==1
        [Y,Z] = meshgrid(params.ts,params.Rs);
        X = ones(size(Y))*t_fixstart;
        s1 = surface(X,Y,Z);
    else
        [X,Z] = meshgrid(params.ts,params.Rs);
        Y = ones(size(X))*t_fixstart;
        s1 = surface(X,Y,Z);
    end
    
    s1.FaceColor = 'k';
    s1.FaceAlpha = 0.5;
    s1.EdgeColor = 'none';
    xyz_this = particle_xyz(i_fixnum,:);
    plot3(xyz_this(:,1),xyz_this(:,2),xyz_this(:,3),'w-','markersize',1,'linewidth',2);
    plot3(xyz_this(1,1),xyz_this(1,2),xyz_this(1,3),'.','Color','w','markersize',20);
    plot3(xyz_this(end,1),xyz_this(end,2),xyz_this(end,3),'.','Color','w','markersize',20);
    plot3(xyz_this(end,1),xyz_this(end,2),xyz_this(end,3),'o','Color','w','markersize',16);
    
    y = 3-y;
    setFigFontSize(28,fh1);
    set(gcf, 'Color', 'w');
    view(-45,20);
end

% Show 2D surface
y = 1;
fh2 = figure('units','normalized','outerposition',[0 0 1 0.6]);
for fn = fixnumunique
    i_fixnum = find(fixnum==fn);
    % End of this fixation overlabs with beginning of next fixation.
    % Take the first timestep of the next fixation
    if i_fixnum(end) < length(fixnum), i_fixnum = [i_fixnum,i_fixnum(end)+1]; end
    
    if fn==1, t_fixstart = 0;
    else, t_fixstart = tItem_all(i_fixnum(1),3-y);
    end
    
    t1 = tItem_all(i_fixnum(1),1);
    t2 = tItem_all(i_fixnum(1),2);
    xyz_this = particle_xyz(i_fixnum,:);
    
    subplot(1,length(fixnumunique),fn);
    if y==1
        opsurface = squeeze(optpolicy.(sprintf('p%d',y))(:,isequal_aij(params.ts,t2,2),:))';
        X = particle_xyz(i_fixnum,2);
    else
        opsurface = squeeze(optpolicy.(sprintf('p%d',y))(isequal_aij(params.ts,t1,2),:,:))';
        X = particle_xyz(i_fixnum,1);
    end
    Z = particle_xyz(i_fixnum,3);

    imagesc(params.ts,params.Rs,opsurface); hold on; axis xy; policyColorSetting;
    plot(X,Z,'w-','markersize',1,'linewidth',0.5);
    plot(X(1),Z(1),'.','Color','w','markersize',20);
    plot(X(end),Z(end),'.','Color','w','markersize',20);
    plot(X(end),Z(end),'o','Color','w','markersize',16);
    set(gca,'ylim',zlimit); pbaspect([1,1,1]); ylabel('Delta(t)');
    if y==1
        set(gca,'XDir','reverse');
        xlabel('t1[s]');
    else
        xlabel('t2[s]');
    end

    y = 3-y;
end





%% Figure 3 & 4B,C
% Behavioral plots for human behavior & optimal policy simulations
%%%%%%%%%%%
% whichTask: 1: aVBDM, 2: aPDM --> Set to 2 to generate Figure 3 - Figure supplement 2
whichTask = 1;
% 0: Human data, 1: Optimal policy simulation
whichModel = 1;
% Plotting options
plotOptions = struct;
plotOptions.saveFigures = 0;

% Other params
showOptPolicy = 0;  % Show optimal policy space
showfixationdurations = 1;  % Show mean fixation duration of the first three fixations
showBehavPlots = 1;  % Show psychometric/chronometric curves
showFixBiasEffects = 1;  % Show fixation bias effects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load behavioral data
taskstr = 'aVBDM';
load(fullfile(data_dir,sprintf('datastruct_%s.mat',taskstr)));
% Load params
params = getParams(whichTask);

% Generate behavior using model
numsubs = length(dstruct_real);
if whichTask==1
    z_reps = 20;  % Number of iterations of trials
    z_all = repmat(permn(0:7,2),z_reps,1);
elseif whichTask==2
    z_reps = 50;  % Number of iterations of trials
    z_all = repmat(permn(0:4,2),z_reps,1);
end
dstruct = getFilledDstruct(dstruct_real,numsubs,z_all);

every10Prc = floor(prctile(1:length(dstruct),10:10:100));
if whichModel == 0  % Human data
    dstruct = dstruct_real;
elseif whichModel == 1  % Optimal model    
    % Load optimal policy information
    if whichTask==1, optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir);
    elseif whichTask==2, optpolicy = getOptPolicy_aPDM_R(params,optpolicy_dir);
    end
    if showOptPolicy == 1
        plotOptions_optPolicy = struct;
        plotOptions_optPolicy.whichTask = 1;
        plotOptions_optPolicy.policySlices = 1;
        plotOptions_optPolicy.policy3DSurf = 1;
        visualizeOptPolicy(optpolicy,params,plotOptions_optPolicy);
    end
    
    % Simulate behavior
    fprintf('Simulating behavior: ');
    for s = 1:length(dstruct)
        if any(s==every10Prc), fprintf('%.0f%% ',find(s==every10Prc)*10); end
        for ti = 1:length(dstruct(s).trialnum)
            [dstruct(s).choice(ti),dstruct(s).rt(ti),dstruct(s).tItem(ti,:),~,dstruct(s).fixdur{ti},dstruct(s).fixitem{ti}]...
                = runTrial_aVBDM(dstruct(s).itemval(ti,:),optpolicy,params);
        end
    end
end

% Plot behavioral plots
if showBehavPlots==1
    bout_sim = getbehavoutput(dstruct);
    stats_behav = makeBehavPlots(bout_sim,plotOptions);
end

% Show fixation bias effects
if showFixBiasEffects==1
    nbin_fn = 4;
    nbin_valsum = 5;
    useFixationBins = 1;
    stats_fixbiaseffects = fixbiaseffects(dstruct,nbin_fn,nbin_valsum,useFixationBins,plotOptions);
end

% Show distribution of fixation durations
if showfixationdurations==1
    numfixshow = 3;
    fixdur_all = cell(1,numfixshow);
    for s = 1:length(dstruct)
        for fn = 1:numfixshow
            fixdur_sub_all = [];
            for ti = 1:length(dstruct(s).fixdur)
                if length(dstruct(s).fixdur{ti}) >= fn
                    fixdur_sub_all = cat(1,fixdur_sub_all,dstruct(s).fixdur{ti}(fn));
                end
            end
            fixdur_sub(s,fn) = mean(fixdur_sub_all);
        end
    end
    
    [mean_fixdur,se_fixdur] = getMeanAndSE(fixdur_sub);
    fh = figure; hold on;
    errorbar(1:numfixshow,mean_fixdur,se_fixdur,'k.','markersize',30);
    set(gca,'xlim',[0,numfixshow+1],'xtick',1:numfixshow);
    ylabel('Fixation duration (s)'); xlabel('Fixation number'); pbaspect([1,1,1]);
end




%% Figure 4D
% Compare mean reward between the aDDM and optimal model

% Load behavioral data
whichTask = 1;  % aVBDM
taskstr = 'aVBDM';
load(fullfile(data_dir,sprintf('datastruct_%s.mat',taskstr)));
% Load params
params = getParams(whichTask);
dt = params.dt;

% Set aDDM parameters
% 1 = Parameters used in original paper (Krajbich et al., 2010)
% 2 = Set parameters to match signal-to-noise ratio with Bayesian model (compare performance)
aDDM_param_type = 2;
if aDDM_param_type==1
    d = 0.2;
    sig2_origModel = 0.02^2;
    aGamma_k = 0.3;
    k = d/dt;
    sig2_k = sig2_origModel/dt;
else
    k = 1;
    aGamma_k = sqrt(params.aGamma);
    sig2_k = 2*params.sig2;
end

% Test aDDM performance for different decision boundaries
decbound_all = 1:1:20;

% Number of subjects & trials
N = length(dstruct_real);
z_trials = 200;
% Get empirical distribution of fixation behavior for each difficulty level
fixdist = getEmpiricalFixationDist(dstruct_real);
maxdectime = 8;

% Load optimal policy information
optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir);

% Store mean accuracy, rt, and reward for both models
acc_addm = nan(N,length(decbound_all));
rt_addm = nan(N,length(decbound_all));
reward_addm = nan(N,length(decbound_all));
acc_opt = nan(N,1);
rt_opt = nan(N,1);
reward_opt = nan(N,1);
every10Prc = floor(prctile(1:N,10:10:100));
for s = 1:N
    if any(s==every10Prc), fprintf('%.0f%% ',find(s==every10Prc)*10); end
    % Generate trials - draw from the prior distribution from the optimal model
    z_all = normrnd(0,sqrt(params.sig2_z),[z_trials,2]);
    dstruct_opt = getFilledDstruct(dstruct_real,1,z_all);
    % Simulate behavior for different boundary conditions
    acc_all = nan(size(z_all,1),length(decbound_all));
    rt_all = nan(size(z_all,1),length(decbound_all));
    reward_all = nan(size(z_all,1),length(decbound_all));
    for i_db = 1:length(decbound_all)
        decbound = decbound_all(i_db);
        dstruct = getFilledDstruct(dstruct_real,1,z_all);
        for ti = 1:length(dstruct.trialnum)
            % Simulate trials for optimal model
            if i_db==1
                [dstruct_opt.choice(ti),dstruct_opt.rt(ti),dstruct_opt.tItem(ti,:),~,dstruct_opt.fixdur{ti},dstruct_opt.fixitem{ti},~]...
                    = runTrial_aVBDM(dstruct_opt.itemval(ti,:),optpolicy,params);
            end
            % Get empirical fixation distribution for the aDDM
            allvaldiff = abs(dstruct.itemval(ti,1)-dstruct.itemval(ti,2));
            if allvaldiff > 8, allvaldiff = 8; end
            fixdist_first = fixdist.all_first{fixdist.valdiff==round(allvaldiff)};
            fixdist_mid = fixdist.all_mid{fixdist.valdiff==round(allvaldiff)};
            
            % first fixation
            if rand <= 0.5, y = 1;
            else, y=2;
            end
            % sequence of fixations
            yseq_toMaxTime = nan(1,maxdectime/dt);
            ni = 1;
            while sum(~isnan(yseq_toMaxTime)) < maxdectime/dt
                if ni==1, itemfixdur = fixdist_first(randperm(length(fixdist_first),1));
                else
                    itemfixdur = fixdist_mid(randperm(length(fixdist_mid),1));
                end
                itemfixN = round(itemfixdur/dt);
                yseq_toMaxTime(ni:ni+itemfixN-1) = y;
                ni = ni+itemfixN;
                y = 3-y;
            end
            yseq_toMaxTime = yseq_toMaxTime(1:maxdectime/dt);
            % Re-write choice
            [dstruct.choice(ti),dstruct.rt(ti),dstruct.fixitem{ti},dstruct.fixdur{ti},dstruct.tItem(ti,:),~] = run_aDDM(dstruct.itemval(ti,:),aGamma_k,dt,k,sig2_k,decbound,yseq_toMaxTime);
        end
        perf = getPerfFromDstruct(dstruct,params);
        acc_all(:,i_db) = perf.acc;
        rt_all(:,i_db) = perf.rt;
        reward_all(:,i_db) = perf.reward;
    end
    acc_addm(s,:) = mean(acc_all,1);
    rt_addm(s,:) = mean(rt_all,1);
    reward_addm(s,:) = mean(reward_all,1);
    % Get rewards earned from optimal model
    perf = getPerfFromDstruct(dstruct_opt,params);
    reward_opt(s) = mean(perf.reward);
    rt_opt(s) = mean(perf.rt);
    acc_opt(s) = mean(perf.acc);
end
fprintf('\n');

% aDDM
[mean_addm_rew,se_addm_rew] = getMeanAndSE(reward_addm);
[mean_addm_acc,se_addm_acc] = getMeanAndSE(acc_addm);
[mean_addm_rt,se_addm_rt] = getMeanAndSE(rt_addm);
% Optimal model
[mean_opt_rew,se_opt_rew] = getMeanAndSE(reward_opt);
[mean_opt_acc,se_opt_acc] = getMeanAndSE(acc_opt);
[mean_opt_rt,se_opt_rt] = getMeanAndSE(rt_opt);

markSize = 32;
% Mean single-trial reward
plots = [];
figure; hold on;
for c = 1:length(mean_addm_rew), plots(1) = errorbar(decbound_all(c),mean_addm_rew(c),se_addm_rew(c),'b.','markersize',markSize); end
% Plot line for optimal policy performance
plots(2) = plot(get(gca,'xlim'),[mean_opt_rew,mean_opt_rew],'r-');
plot(get(gca,'xlim'),[mean_opt_rew-se_opt_rew,mean_opt_rew-se_opt_rew],'r--');
plot(get(gca,'xlim'),[mean_opt_rew+se_opt_rew,mean_opt_rew+se_opt_rew],'r--');
xlabel('Decision boundary height'); ylabel('Mean reward');
legend(plots,{'aDDM','Optimal model'});

% Stats - compare peak reward of aDDM vs optimal perf
perf_addm = reward_addm(:,mean_addm_rew==max(mean_addm_rew));
% perf_addm = reward_addm(:,decbound_all==1);
perf_opt = reward_opt;
[~,p,ci,stats] = ttest2(perf_opt,perf_addm);





%% Figure 4E
% Flexible attention bottleneck

% Load behavioral data
load(fullfile(data_dir,'datastruct_aVBDM.mat'));
% Load params
whichTask = 1;  % VBDM
params = getParams(whichTask);
showOptPolicy = 0;

% Iterate through different values of kappa (attentional modulation of evidence reliability)
kappa_all = [params.aGamma,0.1:0.1:0.9];

% Trials to simulate
N = length(dstruct_real);
z_trials = 5000;
acc_all = nan(N,length(kappa_all));
rew_all = nan(N,length(kappa_all));

% Keep the total Fisher information constant
params_base = params;
sig2_base = params_base.sig2;
aGamma_base = params_base.aGamma;
fisherI = (1+aGamma_base)/sig2_base;
sig2_tot = 1/fisherI;

fprintf('Simulating behavior for all k''s:` ');
every10Prc = floor(prctile(1:N,10:10:100));
for s = 1:N
    if any(s==every10Prc), fprintf('%.0f%% ',find(s==every10Prc)*10); end
    acc_sub = nan(z_trials,length(kappa_all));
    rew_sub = nan(z_trials,length(kappa_all));
    z_all = normrnd(0,sqrt(params.sig2_z),[z_trials,2]);
    for i_k = 1:length(kappa_all)
        if kappa_all(i_k)==params_base.aGamma  % If base model
            params = params_base;
        else
            k = kappa_all(i_k);
            params.sig2 = sig2_tot/(1-k);
            params.aGamma = (params.sig2-sig2_tot)/sig2_tot;
        end
        
        % 1. Load optimal policy information
        optpolicy = getOptPolicy_aVBDM(params,optpolicy_dir);
        if showOptPolicy == 1
            plotOptions_optPolicy = struct;
            plotOptions_optPolicy.policySlices = 1;
            plotOptions_optPolicy.policy3DSurfSimple = 0;
            plotOptions_optPolicy.whichTask = 1;
            visualizeOptPolicy(optpolicy,params,plotOptions_optPolicy);
        end
        
        % 2. Simulate behavior
        choice = nan(size(z_all,1),1);
        for ti = 1:size(z_all,1)
            [choice(ti),~,~,~,~,~,rew_sub(ti,i_k)] = runTrial_aVBDM(z_all(ti,:),optpolicy,params);
        end
        % Get accuracy & reward
        acc_sub(:,i_k) = double(choice == double(z_all(:,1)-z_all(:,2) < 0) + 1);
    end
    rew_all(s,:) = mean(rew_sub,1);
    acc_all(s,:) = mean(acc_sub,1);
end
fprintf('\n');

% Plot mean reward for each kappa
i_red = kappa_all==params_base.aGamma;
i_black = kappa_all~=params_base.aGamma;
[mean_rew_black,se_rew_black] = getMeanAndSE(rew_all(:,i_black));
[mean_rew_red,se_rew_red] = getMeanAndSE(rew_all(:,i_red));
fh1 = figure; hold on;
errorbar(kappa_all(i_black),mean_rew_black,se_rew_black,'k.','markersize',30);
errorbar(kappa_all(i_red),mean_rew_red,se_rew_red,'R.','markersize',30);
set(gca,'xlim',[-0.1,1],'xtick',kappa_all);
% plot(get(gca,'xlim'),[mean_rew_rand,mean_rew_rand],'k--');
xlabel('Attention gamma'); ylabel('Mean reward');
pbaspect([1,0.6,1]);

% Stats - test if there is a peak effect at kappa = 0.5
x = abs(0.5-kappa_all);
b_all = nan(size(rew_all,1),1);
for ti = 1:size(rew_all,1)
    b = regress(rew_all(ti,:)',[ones(length(x),1),x']);
    b_all(ti) = b(2);
end
stats_peakEffect = ttest_full(b_all);




