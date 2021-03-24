function [dx,mean_post] = getPosterior_aVBDM(n,dx,fixseq,tItem,z,params)
%Get the posterior means for two items under the framework of
%attention-based value-based decision making
% Input:
% n: index of time point
% dx: Nx2 matrix of momentary evidence
% ySeq: sequence of attended items

y = fixseq(n);
mean_post = nan(1,2);

% Use flexible gamma?
if ~isfield(params,'aGamma_n')
    % Sequence of attntion weights per time point
    aGammaSeq_1 = params.aGamma.*(fixseq(1:n)==2) + double(fixseq(1:n)==1);
    aGammaSeq_2 = params.aGamma.*(fixseq(1:n)==1) + double(fixseq(1:n)==2);
    % Evidence accumulation
    dx1 = randn*sqrt( (params.sig2/(params.aGamma^(y-1)))*params.dt ) + z(1)*params.dt;
    dx2 = randn*sqrt( (params.sig2/(params.aGamma^(2-y)))*params.dt ) + z(2)*params.dt;
    dx(n,:) = [dx1,dx2];
    % Mean of the posterior distribution
    X1 = sum(dx(1:n,1).*aGammaSeq_1');
    mean_post(1) = (params.zbar*params.sig2 + X1*params.sig2_z) / (params.sig2 + (tItem(1)+params.aGamma*tItem(2))*params.sig2_z);
    X2 = sum(dx(1:n,2).*aGammaSeq_2');
    mean_post(2) = (params.zbar*params.sig2 + X2*params.sig2_z) / (params.sig2 + (params.aGamma*tItem(1)+tItem(2))*params.sig2_z);
else
    aGamma_n = params.aGamma_n;
    aGamma_a = params.aGamma_a;
    % Sequence of attntion weights per time point
    aGammaSeq_1 = aGamma_a.*(fixseq(1:n)==1) + aGamma_n.*(fixseq(1:n)==2);
    aGammaSeq_2 = aGamma_a.*(fixseq(1:n)==2) + aGamma_n.*(fixseq(1:n)==1);
    % Evidence accumulation
    dx1 = randn*sqrt( (params.sig2/(aGamma_n^(y-1) * aGamma_a^(2-y)))*params.dt ) + z(1)*params.dt;
    dx2 = randn*sqrt( (params.sig2/(aGamma_n^(2-y) * aGamma_a^(y-1)))*params.dt ) + z(2)*params.dt;
    dx(n,:) = [dx1,dx2];
    % Mean of the posterior distribution
    X1 = sum(dx(1:n,1).*aGammaSeq_1');
    mean_post(1) = (params.zbar*params.sig2 + X1*params.sig2_z) / (params.sig2 + (tItem(1)*aGamma_a+tItem(2)*aGamma_n)*params.sig2_z);
    X2 = sum(dx(1:n,2).*aGammaSeq_2');
    mean_post(2) = (params.zbar*params.sig2 + X2*params.sig2_z) / (params.sig2 + (tItem(1)*aGamma_n+tItem(2)*aGamma_a)*params.sig2_z);
end

end

