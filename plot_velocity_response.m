function [fig] = plot_velocity_response(data,part,val)
    clc;
    %
    t=data.t;
    u=data.u;
    v1=data.v1;
    v1_filt=data.v1_filt;
    v2=data.v2;
    v2_filt=data.v2_filt;
    %

    figure_title = ['MTE360 - ' part ': ' val];
    
    % Plot data
    fig  = figure('Name',figure_title , 'Color','k', ...
                  'Position',[100 100 1100 800]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    plot_title = strrep('MTE360 – Section 2-1 Overlay for ', '2-1', part);
    plot_title = [plot_title val];
    title(tl,plot_title,'Color','w','FontWeight','normal');

    darken = @(ax) set(ax, 'Color','k','XColor','w','YColor','w', ...
                            'GridColor',[0.35 0.35 0.35], 'MinorGridColor',[0.25 0.25 0.25]);

    ax1 = nexttile; hold(ax1,'on');
    plot(t, v1, 'Color',[0 1 1], 'LineWidth',0.5);              % V1 x in cyan
    plot(t, v1_filt, '--', 'Color',[1 0 1], 'LineWidth',1.8);   % Filtered in magenta
    grid on; xlabel('t [s]'); ylabel('Velocity [mm/s]'); title(['Velocity of Mass 1 – ' part ': ' val],'Color','w');
    legend({'V_1','V_{1 Filtered}'},'Location','eastoutside','TextColor','w','Box','off');
    darken(ax1);
   
    ax2 = nexttile; hold(ax2,'on');
    plot(t, v2, 'Color',[1 0 1], 'LineWidth',0.5);              % X2 x in magenta
    plot(t, v2_filt, '--', 'Color',[0 1 0], 'LineWidth',1.8);   % Filtered in green
    grid on; xlabel('t [s]'); ylabel('Velocity [mm/s]'); title(['Velocity of Mass 1 – ' part ': ' val],'Color','w');
    legend({'V_2','V_{2 Filtered}'},'Location','eastoutside','TextColor','w','Box','off');
    darken(ax2);

    ax3 = nexttile; hold(ax3,'on');
    plot(t, u, 'Color',[0 1 1], 'LineWidth',1.8);           % Lab u in cyan
    grid on; xlabel('t [s]'); ylabel('Input [V]'); title('Control Signal','Color','w');
    darken(ax3);

    linkaxes([ax1,ax2,ax3],'x');
end
