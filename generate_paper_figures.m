%%% FIGURES for astrocyte manuscript

%% Initialization
p = params;
subjects = {'aLeft', 'aRight', 'rLeft', 'rRight', 'sLeft', 'sRight', 'vRight', 'NEWsRight'};
subject_names = {'1L', '1R', '2L', '2R', '3L', '3R', '4R'};

%% FIGURE: Comparison of profile of delay and low-pass filters
%%% load previous setion

FigureName = '2';

time_new = 0:0.001:10;
tau_d = 1.5;
[~, ind] = min(abs(time_new - tau_d));

norm_fac = [1, 1/(tau_d), 1/(tau_d), 4/(sqrt(pi)*tau_d)];   % to keep area under the curve = 1

dirac_max = 0.7;
H_1 = dirac_max*(time_new==time_new(ind))*norm_fac(1);
H_2 = (time_new>=0).*(time_new/tau_d).*exp(-time_new/tau_d)*norm_fac(2);
H_3 = (time_new>=0).*(time_new/tau_d).*exp(-time_new.^2/(2*tau_d^2))*norm_fac(3);
H_4 = (time_new>=0).*(time_new.^2/tau_d^2).*exp(-time_new.^2/(tau_d^2))*norm_fac(4);


fig = figure();
H1 = plot(time_new, H_1, 'k-', 'LineWidth', 2);
hold on;
H2 = plot(time_new, H_2, 'k--', 'LineWidth', 2);
H3 = plot(time_new, H_3, 'k:', 'LineWidth', 2);
H4 = plot(time_new, H_4, 'k-.', 'LineWidth', 2);
hold off;

xlabel('time, $t$ (s)','fontsize',17,'interpreter', 'latex')
ylabel('$h(t)$','fontsize',17,'interpreter', 'latex')
set(gca, 'FontSize', 15)
set(gca, 'xtick',0,'xTickLabel',{'0'},...
    'xLim',[-0.5,10],'ytick',0,'yLim',[0,dirac_max], 'yticklabel',{'0'},'FontSize', 15, 'box', 'on')

text(tau_d-0.1, -0.035, '$\tau_d$', 'fontsize', 20, 'interpreter', 'latex')
% annotation('textbox', [0.25,0.095,0,0], 'String', ...
%         '$\tau_d$', 'fontsize', 15, 'linestyle', 'none', 'interpreter', 'latex')
    
leg = legend(gca, [H1, H2, H3, H4], ...
            ' $h_1(t) = \delta(t - \tau_d)$', ...
            ' $h_2(t) = \frac{t}{\tau_d}e^{-t/\tau_d}$',...
            ' $h_3(t) = \frac{t}{\tau_d}e^{-t^2/2\tau_d^2}$',...
            ' $h_4(t) = \frac{t^2}{\tau_d^2}e^{-t^2/\tau_d^2}$');
set(leg, 'FontSize', 18, 'Orientation','Vertical',...
    'box','off','position',[0.65,0.65,0.01,0.01], 'interpreter', 'latex')

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])


%% Figure: sample experimental response
FigureName = '3';

data = DataExtraction_AstrocyticDelay('aLeft', 5); 

fig = figure('Position', [200, 200, 600, 300]);
cmap1 = colormap('bone');
% cmap1 = colormap('gray');
cmap2 = colormap('hot');
cmap3 = flipud(colormap('gray'));
cmap_new = [cmap1(size(cmap1,1)/2+1:end,:); flipud(cmap2(size(cmap2,1)/2+1:end,:)); cmap3];

for abc = 1
    H = subplot(1,2,1, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    axes('Position', Hpos+[-0.06,0,0.04,0]);

    [~, h2] = contourf(data.experiment.space, data.experiment.time, ...
        data.experiment.BOLDsmooth_all.'-0.5*ones(size(data.experiment.BOLDsmooth_all.')), 12);
    colormap(cmap_new)
    cbar2 = colorbar;
    caxis([-1,1])
    set(gca, 'Fontsize', 12)
    xlabel('$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
    set(cbar2, 'ylim',[-1,0], 'YTick', -1:0.25:0, ....
               'YTickLabel', {-0.5, -0.25, 0, 0.25, 0.5}, 'FontSize', 12)
    ylabel(cbar2, 'BOLD response, $Y(x, t)$','fontsize',15, 'interpreter', 'latex')
    % set(h2, 'LineWidth',0.5)
    set(gca, 'xtick', [-5,-2.5,0,2.5,5], 'xLim',[-5,5],...
             'ytick', 0:5:20,'yLim',[0,20],...
             'FontSize', 12, 'box', 'on')
%     annotation('textbox', get(gca, 'Position')+[-0.09,0.02,0,0], 'String', '(a)', 'fontsize', 18, ...
%             'fontweight', 'b', 'linestyle', 'none')


    H = subplot(1,2,2, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    HAx = axes('Position', Hpos+[0.04,0,0.03,0]);

    plot(data.experiment.time, data.experiment.BOLDraw_center, 'k.', 'MarkerSize',8)
    hold on;
    plot(data.experiment.time, data.experiment.BOLDsmooth_center, 'k-','LineWidth',2)
    hold off;
    xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
    ylabel({'Center response, $Y(x=0, t)$'}, 'fontsize', 15, 'interpreter', 'latex')
    leg = legend('raw', 'smooth');
    set(leg, 'FontSize', 12,'Location','Northeast', 'box', 'off')%, 'Orientation','Vertical','position',[0.905,0.48,0.01,0.10])
    set(gca, 'xtick', 0:5:20, 'xLim',[0,20],...
             'ytick', -0.6:0.2:1, 'yLim', [-0.7,0.8], 'FontSize', 12, 'box', 'on')
%     annotation('textbox', get(gca, 'Position')+[-0.1,0.07,0,0], 'String', '(b)', 'fontsize', 18, ...
%             'fontweight', 'b', 'linestyle', 'none')
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])


%% Figure: Estimated v_b and Gamma across subjects
FigureName = '4';

n_subjects = 7;
relevant_subjects = [1:5, 7];
load singleSourceData/editedFileParams.mat

fig = figure();
h = ploterr(Vmat(relevant_subjects), -Gmat(relevant_subjects),...
        0.3709*Vmat(relevant_subjects)-0.6321, 0.125*(-Gmat(relevant_subjects))+0.1118, ...
         'ko');
set(h, 'MarkerSize', 12, 'MarkerFaceColor', 'k')
% plot(Vmat(1:n_subjects), -Gmat(1:n_subjects), 'ko', 'MarkerFaceColor','k','MarkerSize',12)
set(gca, 'FontSize', 12)
xlabel('wave propagation speed, $\nu_\beta$ (mm s$^{-1}$)','fontsize',15,'interpreter', 'latex')
ylabel('wave damping rate, $\Gamma$ (s$^{-1}$)','fontsize',15,'interpreter', 'latex')
set(gca, 'xtick',0:8,'xLim',[0,8],...
    'ytick',0:0.2:1.6,'yLim',[0,1.6], 'FontSize', 12)
% for i=1:7
%     if Vmat(i)<10
%         text(Vmat(i)+0.4, -Gmat(i), subject_names{i},'FontSize', 15)
%     else
%         text(Vmat(i)-1.2, -Gmat(i), subject_names{i},'FontSize', 15)
%     end
% %     if i<=3
% %         text(Vmat(i)*1.2, -Gmat(i), subject_names{i}) 
% %     elseif i==4 || i==5 || i==7
% %         text(Vmat(i)*1.1, -Gmat(i), subject_names{i}) 
% %     end
% end

text(Vmat(1)+0.2, -Gmat(1)+0.03, subject_names{1},'FontSize', 15)
text(Vmat(2)-0.5, -Gmat(2), subject_names{2},'FontSize', 15)
text(Vmat(3)+0.15, -Gmat(3)+0.06, subject_names{3},'FontSize', 15)
text(Vmat(4)+0.1, -Gmat(4)-0.07, subject_names{4},'FontSize', 15)
text(Vmat(5)-0.35, -Gmat(5)-0.05, subject_names{5},'FontSize', 15)
% text(Vmat(6)+0.5, -Gmat(6)+0.04, subject_names{6},'FontSize', 15)
text(Vmat(7)+0.1, -Gmat(7)-0.05, subject_names{7},'FontSize', 15)
        
set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])


%% Figure: Spatial, temporal and combined stimulus + sample experimental response
FigureName = '5';

x = -5e-3:1e-4:5e-3;
t = 0:0.01:20;
[tt, xx] = meshgrid(t, x);          % mesh matrix of space and time

sigma_x = 1e-3;             % standard deviation of Gaussian spread
zeta_x = exp(-x.^2/(2*sigma_x^2));
stim_on_time = 8;           % how long the stimulus is on
stim_off_time = 12;      % how long the stimulus is off
zeta_t = (t > 0).*(t <= stim_on_time);

% Transform zeta_t_conv into matrix 
zeta_t_mat = repmat(zeta_t, length(x), 1);
zeta_x_mat = repmat(zeta_x.', 1, length(t));

% Get final stimulus z(x,t)
zeta_xt = zeta_x_mat.*zeta_t_mat;

fig = figure('Position', [200, 200, 600, 300]);
% fig = figure();
cmap1 = colormap('bone');
% cmap1 = colormap('gray');
cmap2 = colormap('hot');
cmap3 = flipud(colormap('gray'));
cmap_new = [cmap1(size(cmap1,1)/2+1:end,:); flipud(cmap2(size(cmap2,1)/2+1:end,:)); cmap3];

for abc = 1
    H = subplot(1,3,1, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    axes('Position', Hpos+[-0.05,0.01,-0.01,0]);

    plot(x/1e-3, zeta_x, 'k-', 'LineWidth', 2)
   
    xlabel('$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('spatial profile', 'fontsize', 15, 'interpreter', 'latex')
    set(gca, 'xtick', [-5,-2.5,0,2.5,5], 'xLim',[-5,5],...
             'ytick', 0:0.5:1, 'yLim', [-0.001, 1.2], 'FontSize', 12, 'box', 'on')
    annotation('textbox', get(gca, 'Position')+[0,0.01,0,0], 'String', '(a)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')
    

    H = subplot(1,3,2, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    axes('Position', Hpos+[-0.03,0.01,-0.01,0]);

    plot(t, zeta_t, 'k-', 'LineWidth', 2)
    xlabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('temporal profile, $B(t)$', 'fontsize', 15, 'interpreter', 'latex')
    set(gca, 'xtick', [0.01,8,20], 'xticklabel', {0,8,20},'xLim',[0.01, 20],...
             'ytick', 0:0.5:1, 'yLim', [-0.001, 1.2], 'FontSize', 12, 'box', 'on')
    annotation('textbox', get(gca, 'Position')+[0,0.01,0,0], 'String', '(b)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')


    H = subplot(1,3,3, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    axes('Position', Hpos+[-0.02,0.01,0.07,0]);

    [~,h1] = contourf(x/1e-3, t, zeta_xt.', 10);
    colormap(cmap_new)
    cbar1 = colorbar;
    caxis([-1,1])
    set(gca, 'Fontsize', 12)
    xlabel('$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
    set(cbar1, 'ylim',[0.02,1], 'YTick', [0.02,0.2:0.2:1], ....
               'YTickLabel', {0, 0.2, 0.4, 0.6, 0.8, 1.0}, 'FontSize', 12)
    ylabel(cbar1, 'neural activity, $\phi(x, t)$','fontsize',15, 'interpreter', 'latex')
    set(h1, 'LineWidth',0.2)
    set(gca, 'xtick', [-5,-2.5,0,2.5,5], 'xLim',[-5,5],...
             'ytick', [0.01,8,20], 'yticklabel', {0,8,20},'yLim',[0,20],...
             'FontSize', 12)
    annotation('textbox', get(gca, 'Position')+[0,0.01,0,0], 'String', '(c)', 'fontsize', 18, ...
            'fontweight', 'b', 'linestyle', 'none')
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])

%% Figure: Comparing RAW BOLD signal and models at x=0
% 'model A: old parameters w/o astrocytic delay'
% 'model B: old parameters with astrocytic delay'
% 'model C: new parameters w/o astrocytic delay'
% 'model D: new parameters with astrocytic delay'
FigureName = '6';

letters = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'};
fig = figure('Position', [200, 200, 550, 500]);
for i=[1:5,7]%length(subjects)
    
    data = DataExtraction_AstrocyticDelay(subjects{i}, 5);
    
    H = subplot(4,2,i, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    if mod(i,2)~=0 
        HAx = axes('Position', Hpos+[0, -0.01, 0.05, 0.04]);
    else
        HAx = axes('Position', Hpos+[-0.02, -0.01, 0.05, 0.04]);
    end

    Hexp = plot(data.experiment.time, data.experiment.BOLDraw_center/max(data.experiment.BOLDraw_center(:)), 'ko-', 'MarkerSize',3, 'MarkerFaceColor', 'k');
    hold on;
    H1 = plot(data.experiment.time, data.prediction.models{1}.BOLD_center/max(data.experiment.BOLDraw_center(:)), 'rs-', 'MarkerSize',3, 'MarkerFaceColor', 'r');
    H2 = plot(data.experiment.time, data.prediction.models{2}.BOLD_center/max(data.experiment.BOLDraw_center(:)), 'bd-', 'MarkerSize',3, 'MarkerFaceColor', 'b');
    H3 = plot(data.experiment.time, data.prediction.models{3}.BOLD_center/max(data.experiment.BOLDraw_center(:)), 'm^-', 'MarkerSize',3, 'MarkerFaceColor', 'm');
    H4 = plot(data.experiment.time, data.prediction.models{4}.BOLD_center/max(data.experiment.BOLDraw_center(:)), 'cv-', 'MarkerSize',3, 'MarkerFaceColor', 'c');
    hold off;
    
    set(gca, 'xtick', 0:2:20, 'xLim', [0, 20],...
                 'ytick', -1:0.5:1, 'yLim', [-1.3, 1.2], 'FontSize', 12)
    annotation('textbox', get(gca, 'Position')+[0.325,0.002,0,0], 'String', subject_names{i}, 'fontsize', 18, ...
        'fontweight', 'b', 'linestyle', 'none')
    
    if mod(i,2)==0
        set(gca, 'yticklabel',{})
    end
    
    if i==1
        leg = legend(HAx, [Hexp, H1, H2, H3, H4], 'experiment', ...
            'ORIG\_NOAST', 'ORIG\_AST', 'NEW\_NOAST', 'NEW\_AST');
        set(leg, 'FontSize', 12, ...
            'Orientation','Vertical','position',[0.73,0.25,0.01,0.10], 'box', 'off')
        set(gca, 'xticklabel',{})
    elseif i==2 || i==3  
        set(gca, 'xticklabel',{})
    elseif i==4
        xlabel('time (s)', 'FontSize', 15);
    elseif i==5
        set(gca, 'xticklabel',{})
%         xlabel('time (s)', 'FontSize', 15);
        set(get(gca,'ylabel'),'Position',[-3,1.2],'String','BOLD response','fontsize',15)
    elseif i==7
        xlabel('time (s)', 'FontSize', 15);
    end    
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])

%% Figure: Comparing contour BOLD signal and models 
% 'model A: old parameters w/o astrocytic delay'
% 'model B: old parameters with astrocytic delay'
% 'model C: new parameters w/o astrocytic delay'
% 'model D: new parameters with astrocytic delay'
FigureName = '7';

% letters = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'};
texts = {'experiment', 'ORIG\_NOAST', 'ORIG\_AST', 'NEW\_NOAST', 'NEW\_AST'};
locations = [0.002, 0.013, 0.013, 0.013, 0.013];

cmap1 = colormap('bone');
cmap2 = colormap('hot');
cmap_new = [cmap1(24+1:end,:); flipud(cmap2(24+1:end,:))];

fig = figure('Position', [200, 200, 700, 650]);
for i=[1:5,7]
    data = DataExtraction_AstrocyticDelay(subjects{i}, 5);
    
    tstart_ind = find(data.prediction.time==0);
    tend_ind = find(data.prediction.time==20);
    xstart_ind = find(data.prediction.space==-5*1e-3);
    xend_ind = find(data.prediction.space==5*1e-3);
    
    for j=1:5
        if i==7
            H = subplot(6,5,j+(i-1-1)*5, 'Parent', fig);
        else
            H = subplot(6,5,j+(i-1)*5, 'Parent', fig);
        end
        Hpos = get(H, 'Position');
        delete(H)
        axes('Position', Hpos+[-0.02,0,0.01,0.017]);
        
        
        if j==1
            data_exp = data.experiment.BOLDsmooth_all;
%             [c, h] = contourf(data.experiment.space, data.experiment.time, ...
%                     data_exp.', 12);   % NOT NORMALIZED
            [c, h] = contourf(data.experiment.space, data.experiment.time, ...
                    data_exp.'/max(data_exp(:)), 10);   % NORMALIZED
        else
            data_pred = data.prediction.models{j-1+4}.BOLD_all(xstart_ind:xend_ind,tstart_ind:tend_ind);
            data_pred = data_pred - mean(data_pred(:));     % Correcting fluctuations around mean
%             [c, h] = contourf(data.prediction.space(xstart_ind:xend_ind)/1e-3, ...
%             data.prediction.time(tstart_ind:tend_ind), ...
%             (data_pred.'/max(data_pred(:)))*max(data_exp(:)), 12);  %SCALED TO EXPERIMENT
            [c, h] = contourf(data.prediction.space(xstart_ind:xend_ind)/1e-3, ...
            data.prediction.time(tstart_ind:tend_ind), ...
            (data_pred.'/max(data_pred(:))), 10);    % NORMALIZED
        end
        set(gca, 'Fontsize', 12)
        set(gca, 'xtick', [-5,-2.5,0,2.5,5], 'xLim',[-5,5],...
                 'ytick', 0:5:20,'yLim',[0,20],...
                 'FontSize', 12)  
        colormap(cmap_new)
        caxis([-1,1])
        if i==1
            set(gca, 'xticklabel',{})
            title(texts(j), 'fontsize', 16, 'fontweight', 'b');
%             annotation('textbox', get(gca, 'Position')+[0.045,0.045,0,0], ...
%                 'String', letters{j}, 'fontsize', 18, ...
%             'fontweight', 'b', 'linestyle', 'none')
            cbar = colorbar;
            ylabel(cbar, 'BOLD response','fontsize',15)
            set(cbar, 'Position', [0.91,0.11,0.02,0.83],'FontSize', 12)
        elseif i==7
            xlabel('$x$ (mm)', 'fontsize', 15, 'interpreter', 'latex')
        else
            set(gca, 'xticklabel',{})
        end
        if j==1
            ylabel('$t$ (s)', 'fontsize', 15, 'interpreter', 'latex')
            annotation('textbox', get(gca, 'Position')+[-0.11,-0.03,0,0], 'String', ...
                subject_names{i}, 'fontsize', 18, 'fontweight', 'b', 'linestyle', 'none')
        else
            set(gca, 'yticklabel', {})
        end
        
        % SCALED TO EXPERIMENT                         NORMALIZED
        % aLeft: [-0.43, 0.46] -> [-0.5, 0.5]           -0.93
        % aRight: [-0.57, 0.64] -> [-0.6, 0.7]          -0.89
        % rLeft: [-0.49, 0.68] -> [-0.5, 0.7]           -0.72
        % rRight: [-0.89, 0.64] -> [-0.9, 0.7]          -1.4
        % sLeft: [-0.48, 0.54] -> [-0.5, 0.6]           -0.88
        % sRight: [-0.35, 0.37] -> [-0.4, 0.4]          -0.94
        % vRight: [-0.60, 0.45] -> [-0.6, 0.5]          -1.3
        % NEWsRight: [-0.35, 0.37] -> [-0.4, 0.4]       -0.94
    end
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])

%% OLD FIGURE 7 FOR SPATIAL EXTENT IS REMOVED

%% DATA EXTRACTION: Changes in average parameter estimates for different spatial extent

spatial_extent = [0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 5];

% 4: model A, B, C, D
SummaryEstimates.wf = zeros(4, length(subjects), length(spatial_extent));   
SummaryEstimates.kappa = zeros(4, length(subjects), length(spatial_extent));
SummaryEstimates.tau = zeros(4, length(subjects), length(spatial_extent));
SummaryEstimates.RMSE = zeros(4, length(subjects), length(spatial_extent));
SummaryEstimates.AIC = zeros(4, length(subjects), length(spatial_extent));
for i=1:length(spatial_extent)
    for j=1:8       % subjects  8th subject is NEWsRight
        filename = ['DATA_XT_K_wf=0.1to1_tau=0.1to3_DisregardFirstFiveDataPoints',...
        '/DATA2D_x=',num2str(-spatial_extent(i)),'to',num2str(spatial_extent(i)), '_', ...
                    subjects{j}, '_AstrocyticDelay_DisregardFirstFivePoints.mat'];
        
        load(filename)
        
        for n_models=1:4       % models
            SummaryEstimates.wf(n_models,j,i) = models{n_models}.parameters.wf;
            SummaryEstimates.kappa(n_models,j,i) = models{n_models}.parameters.kappa;
            SummaryEstimates.tau(n_models,j,i) = models{n_models}.parameters.tau;
%             SummaryEstimates.RMSE(n_models,j,i) = models{n_models}.RMSE/sqrt((1+2*(i-1))*75);
            SummaryEstimates.RMSE(n_models,j,i) = models{n_models}.RMSE/sqrt((2*i-1)*75);
            SummaryEstimates.AIC(n_models,j,i) = models{n_models}.AIC;
        end
    end
end
% 
% save('Summary_Estimates_DiffSpatialExtent.mat', 'estimate')

%% FIGURE: Changes in average parameter estimates for different spatial extent
% Run previous section first
FigureName = '8';
spatial_extent = [0, 0.7, 1.4, 2.1, 2.8, 3.5, 4.2, 5];
variables = {'wf', 'kappa', 'tau', 'RMSE', 'AIC'};
original_values = [0.56, 0.65, 0];
ylabels = {'$\omega_f$ (s$^{-1}$)', '$\kappa$ (s$^{-1}$)', ...
    '$\tau_d$ (s)', '$RMSE$'};
ylims = {[0.25, 0.6], [0.15, 0.7], [-0.05, 2.5], [0, 0.7]};
considered_variables = 4;

fig = figure('Position', [200, 200, 500, 550]);
for i=1:considered_variables
    data = eval(cat(2, 'SummaryEstimates.', variables{i}));
    
    if i>3
        error = 1;
    else
        error = sqrt(6);
    end
        
    H = subplot(considered_variables,1,i, 'Parent', fig);
    Hpos = get(H, 'Position');
    delete(H)
    HAx = axes('Position', Hpos+[-0.02, -0.02, 0.05, 0.04]);
    
    
    if i==1 || i==2
        markers = {'ks', 'ko', 'k^-', 'kv-'};
    elseif i==3
        markers = {'ks', 'ko-', 'k^', 'kv-'};
    else
        markers = {'ks-', 'ko-', 'k^-', 'kv-'};
    end
    H1 = errorbar(spatial_extent, squeeze(mean(data(1,[1:5,7],:),2)), ...
        squeeze(std(data(1,[1:5,7],:),0,2))/error, markers{1}, 'MarkerFaceColor','k','MarkerSize',6);
    hold on;
    H2 = errorbar(spatial_extent, squeeze(mean(data(2,[1:5,7],:),2)), ...
        squeeze(std(data(2,[1:5,7],:),0,2))/error, markers{2}, 'Color','r','MarkerFaceColor','r','MarkerSize',6);
    H3 = errorbar(spatial_extent, squeeze(mean(data(3,[1:5,7],:),2)), ...
        squeeze(std(data(3,[1:5,7],:),0,2))/error, markers{3}, 'Color','b','MarkerFaceColor','b','MarkerSize',6);
    H4 = errorbar(spatial_extent, squeeze(mean(data(4,[1:5,7],:),2)), ...
        squeeze(std(data(4,[1:5,7],:),0,2))/error, markers{4}, 'Color','m','MarkerFaceColor','m','MarkerSize',6);
    if i<4
        plot(spatial_extent, original_values(i)*ones(size(spatial_extent)), 'k--', 'LineWidth', 2)
    end
    hold off;
    uistack(H1, 'top')
    set(gca, 'xtick', 0:1:5, 'xLim', [-0.1, 5.1], ...
        'yLim', [ylims{i}(1), ylims{i}(2)], 'FontSize', 12)
    if i==considered_variables
%         xlabel('half of spatial window, $d$ (mm)','fontsize',15, 'interpreter', 'latex')
        xlabel('$d$ (mm)','fontsize',15, 'interpreter', 'latex')
        ylabel(ylabels{i},'fontsize',15,'interpreter','latex')
        leg = legend(HAx, [H1, H2, H3, H4], ...
            'ORIG\_NOAST', 'ORIG\_AST', 'NEW\_NOAST', 'NEW\_AST');
        set(leg, 'FontSize', 12, ...
            'Orientation','Horizontal','position',[0.51,0.97,0.01,0.01], 'box', 'off')
    else
        ylabel(ylabels{i},'fontsize',15,'interpreter','latex')
        set(gca, 'xticklabel', {})
    end
    
    annotation('textbox', [0.005,0.97,0,0], 'String', ...
        '(a)', 'fontsize', 18, 'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.005,0.74,0,0], 'String', ...
        '(b)', 'fontsize', 18, 'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.005,0.53,0,0], 'String', ...
        '(c)', 'fontsize', 18, 'fontweight', 'b', 'linestyle', 'none')
    annotation('textbox', [0.005,0.3,0,0], 'String', ...
        '(d)', 'fontsize', 18, 'fontweight', 'b', 'linestyle', 'none')
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])

%% DATA EXTRACTION: Comparison of data of delay and low-pass filter model at x=0 raw data

% statistics of center only
model_comparison.RMSE_combined = zeros(4, 7, 2);  % 4: delay and low-pass1,2,3, 7: subjects, 2: models with tau as parameter
model_comparison.AIC_combined = zeros(4, 7, 2);  % 4: delay and low-pass1,2,3, 7: subjects, 2: models with tau as parameter

AIC_0 = zeros(1,7);       % AIC of orig_noast for delta AIC

for i=1:length(subjects)
    data = DataExtraction_AstrocyticDelay(subjects{i}, 5); 
    
    model_comparison.subject{i}.experiment_time = data.experiment.time;
    model_comparison.subject{i}.experiment_space = data.experiment.space;   % in mm
    model_comparison.subject{i}.experiment_BOLDraw_center =  data.experiment.BOLDraw_center;
    model_comparison.subject{i}.experiment_BOLDsmooth_all =  data.experiment.BOLDsmooth_all;
    
    model_comparison.subject{i}.prediction_time = data.prediction.time;     
    model_comparison.subject{i}.prediction_space = data.experiment.space;   % in m
    
    AIC_0(i) = data.prediction.models{1}.AIC;
    
    for n_models=[2,4]       % four different models MODEL B, D
        % for delay    CHECK IF models{2} is valid without models{1}
        model_comparison.subject{i}.models{n_models}.prediction_delay_BOLDraw_center =  ...
                                data.prediction.models{n_models}.BOLD_center;
        model_comparison.subject{i}.models{n_models}.prediction_delay_BOLDsmooth_all =  ...
                                data.prediction.models{n_models}.BOLD_all;      %y: distance, x:time
        
        model_comparison.RMSE_combined(1,i,n_models/2) = data.prediction.models{n_models}.RMSE;
        model_comparison.AIC_combined(1,i,n_models/2) = data.prediction.models{n_models}.AIC;
        
        % for low-pass filters
        maxHeight_experiment_center = max(data.experiment.BOLDraw_center);  %% Taking maximum of experiment BOLD
    
        prediction_parameters = [data.prediction.models{n_models}.parameters.wf, ...
                                 data.prediction.models{n_models}.parameters.kappa, ...
                                 data.prediction.models{n_models}.parameters.tau]; 
        
        for lowpass=1:3
            PredictedResponse = HemodynamicResponse1D(prediction_parameters, ...
                            ['low-pass_filter',num2str(lowpass)], p, ...
                            data.v_b*1e-3, data.Gamma);       %% low-pass filter result    
        
            centerResponse_prediction = interp1(PredictedResponse.time, ...
                    real(PredictedResponse.Y_xt(size(PredictedResponse.Y_xt,1)/2+1,:)), ...
                    data.experiment.time);      % interpolate predicted response to experimental time
            centerResponse_prediction = centerResponse_prediction - mean(centerResponse_prediction);    % subtract to measure fluctuations around mean       
            centerResponse_prediction = (centerResponse_prediction/max(centerResponse_prediction))*maxHeight_experiment_center;   %normalize so that experimental and predicted response have the same height
        
            if lowpass==1
                model_comparison.subject{i}.models{n_models}.prediction_lowpass1_BOLDraw_center =  ...
                                        centerResponse_prediction;
                model_comparison.subject{i}.models{n_models}.prediction_lowpass1_BOLDsmooth_all =  ...
                                        real(PredictedResponse.Y_xt);       %y: distance, x:time
            elseif lowpass==2
                model_comparison.subject{i}.models{n_models}.prediction_lowpass2_BOLDraw_center =  ...
                                        centerResponse_prediction;
                model_comparison.subject{i}.models{n_models}.prediction_lowpass2_BOLDsmooth_all =  ...
                                        real(PredictedResponse.Y_xt);       %y: distance, x:time
            elseif lowpass==3
                model_comparison.subject{i}.models{n_models}.prediction_lowpass3_BOLDraw_center =  ...
                                        centerResponse_prediction;
                model_comparison.subject{i}.models{n_models}.prediction_lowpass3_BOLDsmooth_all =  ...
                                        real(PredictedResponse.Y_xt);       %y: distance, x:time
            end

            difference = centerResponse_prediction(5+1:end) - ...
                         model_comparison.subject{i}.experiment_BOLDraw_center(5+1:end);
            n = length(difference);

            model_comparison.RMSE_combined(lowpass+1,i,n_models/2) = sqrt(sum(difference.^2));

            K = n_models-1;
            model_comparison.AIC_combined(lowpass+1,i,n_models/2) = 2*K + n*log((model_comparison.RMSE_combined(lowpass+1,i,n_models/2))^2/n) + ...
                                                            2*K*(K+1)/(n - K - 1);
        end
    end
    display(subjects{i})
end

AIC_0_mat = repmat(AIC_0, 4, 1);
%% FIGURE: Comparison of RMSE of delay and low-pass filter model at x=0 raw data
%%% load previous 4 sections before

FigureName = '9';
n_subjects = 7;
letters = {'(a)', '', '(b)', ''};

fig = figure('Position', [200, 200, 650, 450]);

for i=1:4
    H = subplot(2,2,i);
    Hpos = get(H, 'Position');
    delete(H)

    if i==1
        stats = model_comparison.RMSE_combined(:,[1:5,7],1)/sqrt(75);
        xlabel_name = '';
        ylabel_name = '$RMSE$';
        xtick_label = {};
        yticks = 0:0.1:0.7;
        ylimits = [0, 0.71]; %old limits without sqrt(75) [1, 6];
        ylimits_inset = [0, 0.7];
        adjustment_pos = [-0.01, -0.04, 0.04, 0.04];
        HAx1 = axes('Position', Hpos+adjustment_pos);
    elseif i==2
        stats = model_comparison.AIC_combined(:,[1:5,7],1)-AIC_0_mat(:,[1:5,7],1);
        xlabel_name = '';
        ylabel_name = '$\Delta AIC$';
        xtick_label = {};
        yticks = [-100:40:-20,0,20:40:100];
        ylimits = [-100, 100];
        ylimits_inset = [-100, 10];
        adjustment_pos = [0.02, -0.04, 0.04, 0.04];
        HAx1 = axes('Position', Hpos+adjustment_pos);
    elseif i==3
        stats = model_comparison.RMSE_combined(:,[1:5,7],2)/sqrt(75);
        xlabel_name = 'cases';
        ylabel_name = '$RMSE$';
        xtick_label = subject_names([1:5,7]);
        yticks = 0:0.1:0.7;
        ylimits = [0, 0.71]; %[1, 6];
        ylimits_inset = [0, 0.7];
        adjustment_pos = [-0.01, 0, 0.04, 0.04];
        HAx1 = axes('Position', Hpos+adjustment_pos);
    elseif i==4
        stats = model_comparison.AIC_combined(:,[1:5,7],2)-AIC_0_mat(:,[1:5,7],1);
        xlabel_name = 'cases';
        ylabel_name = '$\Delta AIC$';
        xtick_label = subject_names([1:5,7]);
        yticks = [-100:40:-20,0,20:40:100];
        ylimits = [-100, 100];
        ylimits_inset = [-100, 0];
        adjustment_pos = [0.02, 0, 0.05, 0.05];
        HAx1 = axes('Position', Hpos+adjustment_pos);
    end

    plot(1:6, squeeze(stats(1,:)), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 10)
    hold on;
    plot(1:6, squeeze(stats(2,:)), 'rs', 'MarkerFaceColor','r', 'MarkerSize', 10)
    plot(1:6, squeeze(stats(3,:)), 'b^', 'MarkerFaceColor','b', 'MarkerSize', 10)
    plot(1:6, squeeze(stats(4,:)), 'mv', 'MarkerFaceColor','m', 'MarkerSize', 10)
    if i==2 || i==4
        plot(0:7, zeros(1,8), 'k--', 'LineWidth', 1)
    end
    hold off;
    xlabel(xlabel_name,'fontsize',15)
    ylabel(ylabel_name,'fontsize',15,'interpreter','latex')
    set(gca, 'FontSize', 12)
    set(gca, 'box' ,'off')
    set(gca, 'xtick',1:6,'xTickLabel',xtick_label,...
        'xLim',[0.5,6.5], 'ytick', yticks,'yLim', ylimits, 'FontSize', 12)
    if i==1
        leg = legend('$h_1$', '$h_2$', '$h_3$', '$h_4$');
        set(leg, 'FontSize', 15, 'Orientation','Horizontal',...
            'box','off','position',[0.51,0.98,0.01,0.01], 'interpreter', 'latex')
        text(-0.9, 0.18, 'ORIG\_AST', ...
        'fontsize', 18, 'fontweight', 'b', 'rotation', 90)
    elseif i==3     % old location -0.85, 1,1
        text(-0.9, 0.18, 'NEW\_AST', ...
        'fontsize', 18, 'fontweight', 'b', 'rotation', 90)
    end
%     annotation('textbox', get(gca, 'Position')+[-0.11,0.04,0,0], 'String', letters{i}, 'fontsize', 22, ...
%             'fontweight', 'b', 'linestyle', 'none')

    HAx2 = axes('Position', Hpos+adjustment_pos+[0.06, 0.25, -0.25, -0.25]);
    bar(1, squeeze(mean(stats(1,:),2)), 'FaceColor', 'k')
    hold on;
    bar(2, squeeze(mean(stats(2,:),2)), 'FaceColor', 'r')
    bar(3, squeeze(mean(stats(3,:),2)), 'FaceColor', 'b')
    bar(4, squeeze(mean(stats(4,:),2)), 'FaceColor', 'm')
    errorbar(1:4, squeeze(mean(stats(:,:),2)), ...
        squeeze(std(stats(:,:,1),0,2)),'k.', 'MarkerSize', 1)
    hold off;
%     xlabel('astrocyte functions','fontsize',15)
    ylabel('mean $\pm$ 1sd','fontsize',12,'interpreter','latex')
    set(gca, 'FontSize', 12)
    set(gca, 'box' ,'off')
    set(gca, 'xtick',1:4,'xticklabel', {}, 'xLim',[0.5,4.5], ...
        'ytick',ylimits_inset,'yticklabel',{num2str(ylimits_inset(1)),num2str(ylimits_inset(2))},...
        'ylim',ylimits_inset,'FontSize', 8)
    if i==1 || i==3
        text_loc = 1.55;
    else
        text_loc = 1.25;
    end
    text(0.7, squeeze(mean(stats(1,:),2))+...
        sign(stats(1,1))*text_loc*squeeze(std(stats(1,:),0,2)), '$h_1$', 'fontsize', 11, 'interpreter', 'latex')
    text(1.7, squeeze(mean(stats(2,:),2))+...
        sign(stats(1,1))*text_loc*squeeze(std(stats(2,:),0,2)), '$h_2$', 'fontsize', 11, 'interpreter', 'latex')
    text(2.7, squeeze(mean(stats(3,:),2))+...
        sign(stats(1,1))*text_loc*squeeze(std(stats(3,:),0,2)), '$h_3$', 'fontsize', 11, 'interpreter', 'latex')
    text(3.7, squeeze(mean(stats(4,:),2))+...
        sign(stats(1,1))*text_loc*squeeze(std(stats(4,:),0,2)), '$h_4$', 'fontsize', 11, 'interpreter', 'latex')
end

set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
print('-depsc', ['/Users/jamespang/Desktop/PhD Research/Project2_Astrocyte_Manuscript/Figure',FigureName,'.eps'])

%% TABLE OF PARAMETER VALUES
%%% Load data for each subject
p = params;
subjects = {'aLeft', 'aRight', 'rLeft', 'rRight', 'sLeft', 'sRight', 'vRight', 'NEWsRight'};
subject_names = {'S1L', 'S1R', 'S2L', 'S2R', 'S3L', 'S3R', 'S4R'};

data_aLeft = DataExtraction_AstrocyticDelay('aLeft', 5); 
data_aRight = DataExtraction_AstrocyticDelay('aRight', 5);
data_rLeft = DataExtraction_AstrocyticDelay('rLeft', 5);
data_rRight = DataExtraction_AstrocyticDelay('rRight', 5);    
data_sLeft = DataExtraction_AstrocyticDelay('sLeft', 5);
data_sRight = DataExtraction_AstrocyticDelay('sRight', 5);
data_vRight = DataExtraction_AstrocyticDelay('vRight', 5);
data_NEWsRight = DataExtraction_AstrocyticDelay('NEWsRight', 5);

%%% Table: Parameters with RMSE and AIC 
% RMSE/AIC_all is compared with smooth data, RMSE/AIC_center is with raw data
extent = [-5,5];
model_results = cell(4,1);
relevant_data = [1:5,7];
AIC_all_0 = zeros(1,7);
AIC_center_0 = zeros(1,7);
for i=1:4
    subject = cell(8,1);
    models_wf = zeros(8,2);
    models_kappa = zeros(8,2);
    models_tau = zeros(8,2);
    models_RMSE_all = zeros(8,2);
    models_RMSE_center = zeros(8,2);
    models_AIC_all = zeros(8,2);
    models_AIC_center = zeros(8,2);
    for j=1:7
        data = eval(cat(2, 'data_', subjects{j}));
        sigma2 = diag(svdinv(2*full(data.prediction.models{i}.Jacobian'*data.prediction.models{i}.Jacobian)));
        difference_all = BOLD_LSE_XT([data.prediction.models{i}.parameters.wf, ...
                               data.prediction.models{i}.parameters.kappa, ...
                               data.prediction.models{i}.parameters.tau, 0, 0], ...
                               data.prediction.models{i}.astrocytic_model, p, ...
                               subjects{j}, 5, extent);
        n = length(difference_all);
        K = i - 1;
        
        subject(j) = cellstr(subject_names{j});
        if length(sigma2)==3
            models_wf(j,:) = [data.prediction.models{i}.parameters.wf, sigma2(1)];
            models_kappa(j,:) = [data.prediction.models{i}.parameters.kappa, sigma2(2)];
            models_tau(j,:) = [data.prediction.models{i}.parameters.tau, sigma2(3)];
        elseif length(sigma2)==2
            models_wf(j,:) = [data.prediction.models{i}.parameters.wf, sigma2(1)];
            models_kappa(j,:) = [data.prediction.models{i}.parameters.kappa, sigma2(2)];
            models_tau(j,:) = [data.prediction.models{i}.parameters.tau, 0];
        else
            models_wf(j,:) = [data.prediction.models{i}.parameters.wf, 0];
            models_kappa(j,:) = [data.prediction.models{i}.parameters.kappa, 0];
            models_tau(j,:) = [data.prediction.models{i}.parameters.tau, 0];
        end
        
        Nall = n;
        Ncenter = length(data.experiment.BOLDraw_center)-5;
        models_RMSE_all(j) = sqrt(sum(difference_all.^2))/sqrt(Nall);
        models_RMSE_center(j) = data.prediction.models{i}.RMSE/sqrt(Ncenter);
        models_AIC_all(j) = (2*K + n*log((models_RMSE_all(j))^2*Nall/n) + ...
                               2*K*(K+1)/(n - K - 1))/1000;          % NEED TO MULTIPLY 1000 
        models_AIC_center(j) = data.prediction.models{i}.AIC/1000;     % NEED TO MULTIPLY 1000
        
        if i==1
            AIC_all_0(j) = models_AIC_all(j);
            AIC_center_0(j) = models_AIC_center(j);
        end
        models_AIC_all(j) = models_AIC_all(j) - AIC_all_0(j);
        models_AIC_center(j) = models_AIC_center(j) - AIC_center_0(j);
        
        subject(8) = cellstr('average');
        
        models_wf(8,:) = [mean(models_wf(relevant_data,1)), std(models_wf(relevant_data,1))/sqrt(length(relevant_data))];
        models_kappa(8,:) = [mean(models_kappa(relevant_data,1)), std(models_kappa(relevant_data,1))/sqrt(length(relevant_data))];
        models_tau(8,:) = [mean(models_tau(relevant_data,1)), std(models_tau(relevant_data,1))/sqrt(length(relevant_data))];
        models_RMSE_all(8,:) = [mean(models_RMSE_all(relevant_data)), std(models_RMSE_all(relevant_data))];
        models_RMSE_center(8,:) = [mean(models_RMSE_center(relevant_data)), std(models_RMSE_center(relevant_data))];
        models_AIC_all(8,:) = [mean(models_AIC_all(relevant_data)*1000)/1000, std(models_AIC_all(relevant_data)*1000)/1000];    % NEED TO MULTIPLY 1000 
        models_AIC_center(8,:) = [mean(models_AIC_center(relevant_data)*1000)/1000, std(models_AIC_center(relevant_data)*1000)/1000];    % NEED TO MULTIPLY 1000 
    end

    header = {'w_f','stdev', 'K', 'stdev', 'tau', 'stdev', 'RMSE_all', 'std', 'RMSE_center', 'std', ...
        'AIC_all', 'std', 'AIC_center', 'std'};
    model_results{i} = cat(2,models_wf([1:5,7,8],:), models_kappa([1:5,7,8],:), models_tau([1:5,7,8],:), ...
                        models_RMSE_center([1:5,7,8],:), models_AIC_center([1:5,7,8],:), models_RMSE_all([1:5,7,8],:), models_AIC_all([1:5,7,8],:));
%     table(models_wf([1:5,7,8],:), models_kappa([1:5,7,8],:), models_tau([1:5,7,8],:), ...
%                         models_RMSE_center([1:5,7,8],:), models_AIC_center([1:5,7,8],:), models_RMSE_all([1:5,7,8],:), models_AIC_all([1:5,7,8],:))
end