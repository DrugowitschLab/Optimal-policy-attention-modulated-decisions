function visualizeOptPolicy(optpolicy,params,plotOptions)
%Plot the optimal policy in different ways

yToPlot = 1;

if ~isfield(plotOptions,'saveFigures'), plotOptions.saveFigures = 0; end
if ~isfield(plotOptions,'saveVideos'), plotOptions.saveVideos = 0; end

if plotOptions.whichTask==1
    zlimit = [-5,5];
    zval = params.Rs;
    zStr = '\Delta';
elseif plotOptions.whichTask==2
    zlimit = [0,1];
    zval = params.gs;
    zStr = 'g';
else, error('Specify decision type (aVBDM or aPDM)');
end

% Settings
colors_policy = linspecer(4); colors_y = [1,0,0;0,0,1];
% Video settings
figvidOption.FrameRate=30; figvidOption.Duration=5.5; figvidOption.Periodic=true; figvidOption.LightEveryAngle = true;
vidViews = [0,30;90,30;180,30;270,30;359,30];
% vidViews = [-20,10;-110,10;-190,80;-290,10;-380,10];
saveStrAppend = getOptPolicySaveStrAppend(params);

% Get boundaries for optimal policy
optbound = getOptpolicyBoundary(optpolicy,params);

% Visualize optimal policy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlimit = [min(params.ts),max(params.ts)]; ylimit = [min(params.ts),max(params.ts)];

if isfield(plotOptions,'policySlices') && plotOptions.policySlices==1
    % Show optimal policy across time
    % (For fixed t2, plot policy vs t1 and vice versa)
    n_increments = 1:20:100;
    spIdx = [1,2];
    y = yToPlot;
    fh = figure('units','normalized','outerposition',[0 0 0.4 1]);
    for t_show = 1:2
        for ni = 1:length(n_increments)
            n_other = n_increments(ni);
            if t_show==1, p1_t = squeeze(optpolicy.(sprintf('p%d',y))(:,n_other,:))';
            else, p1_t = squeeze(optpolicy.(sprintf('p%d',y))(n_other,:,:))';
            end
            subplot(length(n_increments),2,spIdx(t_show));
            imagesc(params.ts,zval,p1_t); axis xy; policyColorSetting;
            xlabel(sprintf('t%d (s)',t_show)); ylabel(zStr); title(sprintf('y=%d (t%d = %.01f)',y,3-t_show,params.ts(n_other)));
            spIdx(t_show) = spIdx(t_show)+2;
        end
    end
    if plotOptions.saveFigures == 1, print(fh,fullfile(plotOptions.figSaveDir,sprintf('optPolicySlices1%s.tiff',saveStrAppend)),'-dtiff','-r100'); end
    
    % Plot top-down slices of the decision boundaries
    if plotOptions.whichTask==1
        z_toShow = -3:1.5:3;
    elseif strcmp(plotOptions.type,'aPDM'), z_toShow = [0.3,0.4,0.5,0.6,0.7];
    end
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    spIdx = 1;
    for y = 1:2
        for i = 1:length(z_toShow)
            R_i = zval==z_toShow(i);
            subplot(3,length(z_toShow),spIdx); imagesc(params.ts,params.ts,optpolicy.(sprintf('p%d',y))(:,:,R_i)); hold on; axis xy; policyColorSetting;
            title(sprintf('y = %d (%s = %.1f)',y,zStr,z_toShow(i))); xlabel('t2'); ylabel('t1');
            pbaspect([1,1,1]);
            spIdx = spIdx+1;
        end
    end
    lw = 1;
    for ri = 1:length(z_toShow)
        makePlots = true;
        r = z_toShow(ri);
        plots = nan(2,1);
        subplot(3,length(z_toShow),ri+2*length(z_toShow)); hold on;
        for y=1:2
            z_i_dec = optbound.dec{y}(:,3)==r;
            z_i_accum = optbound.accum{y}(:,3)==r;
            if ~any(z_i_dec) && ~any(z_i_accum), makePlots = false;
            else
                plot(params.ts(optbound.dec{y}(z_i_dec,1)), params.ts(optbound.dec{y}(z_i_dec,2)),'color', colors_y(y,:), 'LineWidth', lw);
                plots(y) = plot(params.ts(optbound.accum{y}(z_i_accum,1)), params.ts(optbound.accum{y}(z_i_accum,2)),'color', colors_y(y,:), 'LineWidth', lw);
            end
            % There may by no switch boundaries
            if ~isempty(optbound.switch{y})
                z_i_switch = optbound.switch{y}(:,3)==r;
                plot(params.ts(optbound.switch{y}(z_i_switch,1)), params.ts(optbound.switch{y}(z_i_switch,2)),'color', colors_y(y,:), 'LineWidth', lw);
            end
        end
        if ri==1 && makePlots, legend(plots,{'y = 1','y = 2'}); end
        title(sprintf('%s = %.1f',zStr,r)); set(gca,'xlim',[0,params.ts(end)],'ylim',[0,params.ts(end)]);
        xlabel('t2'); ylabel('t1'); pbaspect([1,1,1]);
    end
    if plotOptions.saveFigures == 1
%         print(fh,fullfile(plotOptions.figSaveDir,sprintf('optPolicySlices2%s.eps',saveStrAppend)),'-depsc','-r200');
        print(fh,fullfile(plotOptions.figSaveDir,sprintf('optPolicySlices2%s.tiff',saveStrAppend)),'-dtiff','-r100');
    end
end

% Show volumes filled with policy space
if isfield(plotOptions,'policy3DVol') && plotOptions.policy3DVol==1
    % 3d grid settings
    xspacing = xlimit(1):xlimit(end);
    yspacing = ylimit(1):ylimit(end);
    zspacing = zlimit(1):zlimit(end);
    linecolor = 'k';
    linealpha = 0.1;
    
    allAxes3 = [];
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    for pp = 1:4
        alphaMat = double(optpolicy.p1==pp);
        allAxes3(pp) = subplot(1,4,pp);
        vol3d('CData',optpolicy.p1,'XData',params.ts,'YData',params.ts,'ZData',zval,'Alpha',alphaMat);
        colormap(linspecer);
        plot3dgrid(xspacing,yspacing,zspacing,linecolor,linealpha);
        set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit,'xtick',xlimit(1):2:xlimit(end),'ytick',ylimit(1):2:ylimit(end),'ztick',zlimit(1):zlimit(end));
        ylabel('t1'); xlabel('t2');
        if plotOptions.whichTask==1, pbaspect([range(xlimit),range(ylimit),range(zlimit)]); end
        camlight(-30,10);
    end
    Link = linkprop(allAxes3,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    savefilename = sprintf('optPolicyVol%s',saveStrAppend);
    if plotOptions.saveVideos==1, CaptureFigVid(vidViews, fullfile(plotOptions.figSaveDir,savefilename),figvidOption); end
    
    % plot all policies in one plot (excluding choose 1)
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    for pp = 2:4
        alphaMat = double(optpolicy.p1==pp);
        vol3d('CData',optpolicy.p1,'XData',params.ts,'YData',params.ts,'ZData',zval,'Alpha',alphaMat);
        colormap(linspecer);
        plot3dgrid(xspacing,yspacing,zspacing,linecolor,linealpha);
        set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit,'xtick',xlimit(1):2:xlimit(end),'ytick',ylimit(1):2:ylimit(end),'ztick',zlimit(1):zlimit(end));
        ylabel('t1'); xlabel('t2');
        if plotOptions.whichTask==1, pbaspect([range(xlimit),range(ylimit),range(zlimit)]);
        elseif plotOptions.whichTask==1, pbaspect([1,1,1]);
        end
    end
end

% Plot reconstructed surface
if isfield(plotOptions,'policy3DSurf') && plotOptions.policy3DSurf==1
    fh = figure;
    allAxes = [];
    axison = 1;
    deltastr = '\Delta';
    for y = 1:2
        allAxes(y) = subplot(1,2,y); hold on;
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
                %         shading interp; colormap(linspecer);
                grid on;
                if axison==0, axis off; end
                light_handle = light;
                lightangle(light_handle,-30,0);
                material dull;
                set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit);
                ylabel('t1[s]'); xlabel('t2[s]'); zlabel(sprintf('%s(t)',deltastr))
                pbaspect([range(xlimit),range(ylimit),range(zlimit)]);
                
                %         % 3D grid
                %         xspacing = xlimit(1):xlimit(end);
                %         yspacing = ylimit(1):ylimit(end);
                %         zspacing = zlimit(1):zlimit(end);
                %         linecolor = 'k';
                %         linealpha = 0.1;
                %         plot3dgrid(xspacing,yspacing,zspacing,linecolor,linealpha);
            end
        end
        set(gca,'xtick',0:2:6,'ytick',0:2:6,'ztick',-2:2);
        setFigFontSize(26,fh);
        lightangle(light_handle,80,0);
    end
    Link = linkprop(allAxes,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
end


%%%%%%%%%%%%%%%%%%%%%%%% OLDER METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(plotOptions,'policy3DSlices') && plotOptions.policy3DSlices==1
    % Show 3D plot with different slices
    y = yToPlot;
    if plotOptions.whichTask==1
        t_slices = [0,1,2]; z_slices = [0,1];
        z_min = -8; z_max = 8; z_incr = 2;
    elseif plotOptions.whichTask==2
        t_slices = [0,1,2]; z_slices = [0.5,0.7];
        z_min = 0.2; z_max = 0.8; z_incr = 0.2;
    end
    zs_use_i = zval >= z_min & zval <= z_max;
    spIdx = 1:length(t_slices); allAxes1 = []; axn = 1;
    
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    for col = 1:length(t_slices)
        for row = 1:length(z_slices)
            slice_t1 = t_slices(col); slice_t2 = t_slices(col);
            slice_R = z_slices(row);
            allAxes1(axn) = subplot(length(z_slices),length(t_slices),spIdx(col));
            s = slice(params.ts,params.ts,zval(zs_use_i),optpolicy.(sprintf('p%d',y))(:,:,zs_use_i),slice_t1,slice_t2,slice_R); hold on;
            set(s,'edgecolor','none','facealpha',0.8);
            set(gca,'xlim',[0,params.ts(end)],'xtick',0:params.ts(end),'ylim',[0,params.ts(end)],'ytick',0:params.ts(end),'zlim',[z_min,z_max],'ztick',z_min:z_incr:z_max);
            plotMeshLines(s(1),5,5); plotMeshLines(s(2),5,5);
            policyColorSetting; pbaspect([1,1,0.8]);
            xlabel('t2'); ylabel('t1'); zlabel(zStr); title(sprintf('y = %d',y));
            spIdx(col) = spIdx(col)+length(t_slices); axn = axn+1;
        end
    end
    Link = linkprop(allAxes1,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    if plotOptions.saveFigures == 1
        savefig(fh,fullfile(plotOptions.figSaveDir,sprintf('optPolicySlices3D%s.fig',saveStrAppend)));
    end
end

if isfield(plotOptions,'policy3DSurf_old') && plotOptions.policy3DSurf==1
    % Plot 3D surface for switch boundaries and decision boundary for a single
    % item attention
    y = yToPlot; allAxes2 = [];
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    fNames = {'dec','accum','switch'};
    policyColors = linspecer(length(fNames)); axCount = 1;
    for pp = 1:2
        if pp==2
            allAxes2(axCount) = subplot(1,length(fNames)+1,length(fNames)+1); axCount = axCount+1;
            title('All policies');
        end
        for fi = 1:length(fNames)
            if ~isempty(optbound.(fNames{fi}){y}) % a policy might not exist (e.g. no switch boundaries)
                R_i = ~isnan(optbound.(fNames{fi}){y}(:,3));
                % Remove duplicate points
                xx = params.ts(optbound.(fNames{fi}){y}(R_i,1));
                yy = params.ts(optbound.(fNames{fi}){y}(R_i,2));
                zz = optbound.(fNames{fi}){y}(R_i,3);
                if pp==1
                    allAxes2(axCount) = subplot(1,length(fNames)+1,fi); axCount = axCount+1;
                    title(sprintf('y=%d (%s)',y,fNames{fi}));
                end
                hold on;
                scatter3(xx,yy,zz,1,policyColors(fi,:),'.','LineWidth',0.01); view([-42,35]);
                set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit);
                ylabel('t1'); xlabel('t2'); pbaspect([1,1,1]);
            end
        end
    end
    Link = linkprop(allAxes2,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    savefilename = sprintf('optPolicySurf1%s',saveStrAppend);
    if plotOptions.saveFigures == 1, savefig(fh,fullfile(plotOptions.figSaveDir,sprintf('%s.fig',savefilename))); end
    if plotOptions.saveVideos==1, CaptureFigVid(vidViews, fullfile(plotOptions.figSaveDir,savefilename),figvidOption); end

    
    % Plot reconstructed surface
    y = yToPlot; allAxes2 = [];
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    fNames = {'dec','accum','switch'};
    for pp = 1:2
        if pp==2
            allAxes2(axCount) = subplot(1,length(fNames)+1,length(fNames)+1); axCount = axCount+1;
            title('All policies');
        end
        for fi = 1:length(fNames)
            if ~isempty(optbound.(fNames{fi}){y}) % a policy might not exist (e.g. no switch boundaries)
                R_i = ~isnan(optbound.(fNames{fi}){y}(:,3));
                % Remove duplicate points
                xx = params.ts(optbound.(fNames{fi}){y}(R_i,1))';
                yy = params.ts(optbound.(fNames{fi}){y}(R_i,2))';
                zz = optbound.(fNames{fi}){y}(R_i,3);
                if pp==1
                    allAxes2(axCount) = subplot(1,length(fNames)+1,fi); axCount = axCount+1;
                    title(sprintf('y=%d (%s)',y,fNames{fi}));
                end
                hold on;
                
                % Run  program
                trisurface=MyCrustOpen([xx,yy,zz]);
                
                % plot the output triangulation
                trisurf(trisurface,xx,yy,zz,'facecolor','b','edgecolor','b','facealpha',0.4);
                shading interp; colormap(linspecer); grid on;
                lightangle(-45,0); material metal;
                set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit);
                ylabel('t1'); xlabel('t2');
                if plotOptions.whichTask==1, pbaspect([range(xlimit),range(ylimit),range(zlimit)]);
                elseif plotOptions.whichTask==1, pbaspect([1,1,1]);
                end
                
                xspacing = xlimit(1):xlimit(end);
                yspacing = ylimit(1):ylimit(end);
                zspacing = zlimit(1):zlimit(end);
                linecolor = 'k';
                linealpha = 0.1;
                plot3dgrid(xspacing,yspacing,zspacing,linecolor,linealpha);
            end
        end
    end
    Link = linkprop(allAxes2,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    savefilename = sprintf('optPolicySurfRecon1_%s',saveStrAppend);
    if plotOptions.saveFigures == 1, savefig(fh,fullfile(plotOptions.figSaveDir,sprintf('%s.fig',savefilename))); end
    if plotOptions.saveVideos==1, CaptureFigVid(vidViews, fullfile(plotOptions.figSaveDir,savefilename),figvidOption); end
    
    
    % Plot 3D surface of the boundaries in both value functions
    allAxes3 = [];
    fh = figure('units','normalized','outerposition',[0 0 1 1]);
    vfuncolors = ([1,0,0;0,0,1]);
    for fi = 1:length(fNames)
        allAxes3(fi) = subplot(1,length(fNames),fi); hold on;
        plots = [];
        for y = 1:2
            if ~isempty(optbound.(fNames{fi}){y})
                R_i = ~isnan(optbound.(fNames{fi}){y}(:,3));
                xx = params.ts(optbound.(fNames{fi}){y}(R_i,1));
                yy = params.ts(optbound.(fNames{fi}){y}(R_i,2));
                zz = optbound.(fNames{fi}){y}(R_i,3);
                plots(y) = scatter3(xx,yy,zz,1,vfuncolors(y,:),'.','LineWidth',0.01); view([-42,35]);
                title(fNames{fi}); ylabel('t1'); xlabel('t2'); pbaspect([1,1,1]);
                set(gca,'xlim',[min(params.ts),max(params.ts)],'ylim',[min(params.ts),max(params.ts)]);
            end
        end
        legend(plots,{'y=1','y=2'});
    end
    Link = linkprop(allAxes3,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    savefilename = sprintf('optPolicySurf2%s',saveStrAppend);
    if plotOptions.saveFigures == 1, savefig(fh,fullfile(plotOptions.figSaveDir,sprintf('%s.fig',savefilename))); end
    if plotOptions.saveVideos==1, CaptureFigVid(vidViews, fullfile(plotOptions.figSaveDir,savefilename),figvidOption); end
end

if isfield(plotOptions,'policy3DSurfSimple') && plotOptions.policy3DSurfSimple==1
    y = yToPlot;
    fh = figure; hold on; title('All policies');
    fNames = {'dec','accum','switch'};
    for fi = 1:length(fNames)
        if ~isempty(optbound.(fNames{fi}){y}) % a policy might not exist (e.g. no switch boundaries)
            R_i = ~isnan(optbound.(fNames{fi}){y}(:,3));
            % Remove duplicate points
            xx = params.ts(optbound.(fNames{fi}){y}(R_i,1))';
            yy = params.ts(optbound.(fNames{fi}){y}(R_i,2))';
            zz = optbound.(fNames{fi}){y}(R_i,3);
            
            % Run  program
            trisurface=MyCrustOpen([xx,yy,zz]);
            
            % plot the output triangulation
            trisurf(trisurface,xx,yy,zz,'facecolor','b','edgecolor','b','facealpha',0.4);
            shading interp; colormap(linspecer); grid on;
            lightangle(-45,0); material metal;
            set(gca,'xlim',xlimit,'ylim',ylimit,'zlim',zlimit);
            ylabel('t1'); xlabel('t2');
            if plotOptions.whichTask==1, pbaspect([range(xlimit),range(ylimit),range(zlimit)]);
            elseif plotOptions.whichTask==1, pbaspect([1,1,1]);
            end
            
            xspacing = xlimit(1):xlimit(end);
            yspacing = ylimit(1):ylimit(end);
            zspacing = zlimit(1):zlimit(end);
            linecolor = 'k';
            linealpha = 0.1;
            plot3dgrid(xspacing,yspacing,zspacing,linecolor,linealpha);
        end
    end
    savefilename = sprintf('optPolicySurfRecon1_%s',saveStrAppend);
    if plotOptions.saveFigures == 1, savefig(fh,fullfile(plotOptions.figSaveDir,sprintf('%s.fig',savefilename))); end
    if plotOptions.saveVideos==1, CaptureFigVid(vidViews, fullfile(plotOptions.figSaveDir,savefilename),figvidOption); end
end



checksymmetry = 0;
if checksymmetry == 1
    policyEqual = zeros(size(optpolicy.p1,3),1);
    for i_r = 1:size(optpolicy.p1,3)
        p1 = optpolicy.p1(:,:,i_r);
        p2 = optpolicy.p2(:,:,i_r);
        p1(isnan(p1))=3; p2(isnan(p2))=3;
        policyEqual(i_r) = isequal(p1,p2');
    end
end


end

