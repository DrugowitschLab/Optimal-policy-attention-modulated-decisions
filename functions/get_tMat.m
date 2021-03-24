function tMat = get_tMat(Rs,t1,t2,params)
%Get transition matrix T = p(R(t+dt)|R(t),t1,t2,y)
% x-axis: R(t)
% y-axis: R(t+dt)
% z-axis: y (attended item)

dt = params.dt;
sig2_z = params.sig2_z; sig2 = params.sig2;
tMat = nan(length(Rs),length(Rs),2);
tMat_var = nan(2,1);

aGamma = params.aGamma;
sig2_r1_y1 = (sig2*dt) / ( (sig2/sig2_z)+t1+aGamma*t2+dt )^2;
sig2_r2_y1 = (sig2*aGamma*dt) / ( (sig2/sig2_z)+aGamma*t1+t2+aGamma*dt )^2;
sig2_r1_y2 = (sig2*aGamma*dt) / ( (sig2/sig2_z)+t1+aGamma*t2+aGamma*dt )^2;
sig2_r2_y2 = (sig2*dt) / ( (sig2/sig2_z)+aGamma*t1+t2+dt )^2;

tMat_var(1) = (sig2_r1_y1 + sig2_r2_y1)/4;
tMat_var(2) = (sig2_r1_y2 + sig2_r2_y2)/4;

for y = 1:2
    tMat(:,:,y) = normpdf(Rs,Rs',sqrt(tMat_var(y)));
    tMat(:,:,y) = bsxfun(@rdivide, tMat(:,:,y), sum(tMat(:,:,y), 1));
end

end

