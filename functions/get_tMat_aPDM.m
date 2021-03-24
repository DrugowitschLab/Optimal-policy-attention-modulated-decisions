function tMat = get_tMat_aPDM(params,t1,t2)
%Get transition matrix T = p(R(t+dt)|R(t),t1,t2,y)
% x-axis: g
% y-axis: g'
% z-axis: y (attended item)

dt = params.dt;
sig2_z = params.sig2_z; sig2 = params.sig2;
aGamma = params.aGamma;
g = params.gs;

tMat = nan(length(g),length(g),2);
s = sig2/sig2_z;
sig2_r = params.sig2/(s+t1+aGamma*t2) + sig2/(s+aGamma*t1+t2);
for y = 1:2
    sig2_R = ( ( (sig2*(aGamma^(y-1))*dt) / (s+t1+aGamma*t2+dt*aGamma^(y-1))^2 ) + ( (sig2*(aGamma^(2-y))*dt) / (s+aGamma*t1+t2+dt*aGamma^(2-y))^2 ) )/4;
    t1_next = t1+dt*(y==1); t2_next = t2+dt*(y==2);
    sig2_r_next = params.sig2/(s+t1_next+aGamma*t2_next) + sig2/(s+aGamma*t1_next+t2_next);
    tMat(:,:,y) = sqrt(sig2_r_next/sig2_R) * exp( ( ((sig2_R-sig2_r_next)*(norminv(g').^2)) + 2*sqrt(sig2_r_next*sig2_r)*norminv(g')*norminv(g) - sig2_r*(norminv(g).^2)) / (2*sig2_R) );
    tMat(:,:,y) = bsxfun(@rdivide, tMat(:,:,y), sum(tMat(:,:,y), 1));
end


% % % Plot transition matrices for y = 1,2
% % figure;
% % for y = 1:2
% %     subplot(1,2,y);
% %     imagesc(params.gs,params.gs,tMat(:,:,y)); colorbar; colormap(linspecer);
% %     xlabel('g'); ylabel('g''');
% %     pbaspect([1,1,1]);
% % end
% % 
% % 
% % % Sanity check
% % % p(g'|g) should be a martingale
% % t1 = 1; t2 = 1;
% % sig2_z = params.sig2_z; sig2 = params.sig2;
% % aGamma = params.aGamma;
% % g = params.gs;
% % 
% % dt_all = [0.1,0.05,0.01,0.005,0.001];
% % colors = linspecer(length(dt_all));
% % plots = [];
% % figure; hold on;
% % for i_dt = 1:length(dt_all)
% %     dt = dt_all(i_dt);
% %     tMat_this = nan(length(g),length(g));
% %     s = sig2/sig2_z;
% %     sig2_r = params.sig2/(s+t1+aGamma*t2) + sig2/(s+aGamma*t1+t2);
% %     y = 1;
% %     sig2_R = (sig2*(aGamma^(y-1))*dt) / (s+t1+aGamma*t2+dt*aGamma^(y-1))^2 + (sig2*(aGamma^(2-y))*dt) / (s+aGamma*t1+t2+dt*aGamma^(2-y))^2;
% %     t1_next = t1+dt*(y==1); t2_next = t2+dt*(y==2);
% %     sig2_r_next = params.sig2/(s+t1_next+aGamma*t2_next) + sig2/(s+aGamma*t1_next+t2_next);
% %     tMat_this = sqrt(sig2_r_next/sig2_R) * exp( (((sig2_R-sig2_r_next)*norminv(g').^2) + 2*sqrt(sig2_r*sig2_r_next)*norminv(g')*norminv(g) - sig2_r*norminv(g).^2) / (2*sig2_R) );
% %     tMat_this = bsxfun(@rdivide, tMat_this, sum(tMat_this, 1));
% %     
% %     temp = [];
% %     for i = 1:size(tMat_this,2)
% %         temp = cat(1,temp,sum(tMat_this(:,i).*g'));
% %     end
% %     plots(i_dt) = plot(g,g-temp','Color',colors(i_dt,:),'linewidth',2);
% % end
% % legend(plots,sprintfc('%.03f',dt_all));
% % ylabel('Error (E(g''|g)-g)'); xlabel('g');
% % 
% % % Example plots
% % 
% % t1 = 5; t2 = 0.05;
% % sig2_z = params.sig2_z; sig2 = params.sig2;
% % aGamma = params.aGamma;
% % g = params.gs;
% % 
% % dt_all = [0.2,0.05,0.001];
% % colors = linspecer(length(dt_all));
% % plots = [];
% % figure; hold on;
% % for i_dt = 1:length(dt_all)
% %     dt = dt_all(i_dt);
% %     tMat_this = nan(length(g),length(g));
% %     s = sig2/sig2_z;
% %     sig2_r = params.sig2/(s+t1+aGamma*t2) + sig2/(s+aGamma*t1+t2);
% %     y = 1;
% %     sig2_R = (sig2*(aGamma^(y-1))*dt) / (s+t1+aGamma*t2+dt*aGamma^(y-1))^2 + (sig2*(aGamma^(2-y))*dt) / (s+aGamma*t1+t2+dt*aGamma^(2-y))^2;
% %     t1_next = t1+dt*(y==1); t2_next = t2+dt*(y==2);
% %     sig2_r_next = params.sig2/(s+t1_next+aGamma*t2_next) + sig2/(s+aGamma*t1_next+t2_next);
% %     tMat_this = sqrt(sig2_r_next/sig2_R) * exp( (((sig2_R-sig2_r_next)*norminv(g').^2) + 2*sqrt(sig2_r*sig2_r_next)*norminv(g')*norminv(g) - sig2_r*norminv(g).^2) / (2*sig2_R) );
% %     tMat_this = bsxfun(@rdivide, tMat_this, sum(tMat_this, 1));
% %     subplot(1,3,i_dt);
% %     imagesc(params.gs,params.gs,tMat_this); colorbar; pbaspect([1,1,1]); colormap(linspecer);
% %     title(sprintf('dt=%.03f (t1=%.02f,t2=%.02f)',dt,t1,t2));
% %     xlabel('g'); ylabel('g''');
% % end




end

