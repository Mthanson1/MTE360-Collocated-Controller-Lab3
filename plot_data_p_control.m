function [fig] = plot_data_p_control(data,part,kp)
    clc;
    %Load Data
    % Check if data has Scope_DataX1 field
    if isfield(data, 'Scope_DataX1') && isfield(data.Scope_DataX1, 'value')
        t_exp=data.ScopeData_x1.time;
        xr=data.ScopeData_xr.signals.values;
        x1=data.ScopeData_x1.signals.values;
        x2=data.ScopeData_x2.signals.values;
        u=data.ScopeData_u.signals.values;
    elseif isfield(data, 't_exp')
        % Fallback if Scope_DataX1 doesn't exist
        t_exp = data.t_exp;
        xr    = data.xr;
        x1    = data.x1;
        x2    = data.x2;
        u     = data.u;
    else
        error('No suitable time data found in this trial.');
    end
    
    figure_title = ['MTE360 - ' part];
    

    % Plot data
    fig  = figure('Name',figure_title , 'Color','k', ...
                  'Position',[100 100 1100 800]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    plot_title = strrep('MTE360 – Section 2-1 Overlay for K_p= ', '2-1', part);
    plot_title = [plot_title num2str(kp)];
    title(tl,plot_title,'Color','w','FontWeight','normal');

    darken = @(ax) set(ax, 'Color','k','XColor','w','YColor','w', ...
                            'GridColor',[0.35 0.35 0.35], 'MinorGridColor',[0.25 0.25 0.25]);

    ax1 = nexttile; hold(ax1,'on');
    plot(t_exp, xr, '--', 'Color',[0 1 0], 'LineWidth',1.6);   % reference in green
    plot(t_exp, x1, 'Color',[0 1 1], 'LineWidth',1.8);         % X1 x in cyan
    plot(t_exp, x2, 'Color',[1 0 1], 'LineWidth',1.2);         % X2 x in magenta
    grid on; xlabel('t [s]'); ylabel('Position [mm]'); title(['Position – ' part],'Color','w');
    legend({'X_r','X_1','X_2'},'Location','eastoutside','TextColor','w','Box','off');
    darken(ax1);

    %f(1)=subplot(3,1,1);
    %plot(t_exp,xr,t_exp,x1,t_exp,x2)
    %ylabel('Position [mm]')
    %legend('Xr','X1','X2')
    %grid on
   
    ax2 = nexttile; hold(ax2,'on');
    plot(t_exp, xr-x2, 'Color',[0 1 1], 'LineWidth',1.8);     % Lab err in cyan
    grid on; xlabel('t [s]'); ylabel('Error [mm]'); title('Tracking Error X_r → X_2','Color','w');
    darken(ax2);

    %f(2)=subplot(3,1,2);
    %plot(t_exp,xr-x2);
    %ylabel('Error [mm]')
    %grid on;

    ax3 = nexttile; hold(ax3,'on');
    plot(t_exp, u, 'Color',[0 1 1], 'LineWidth',1.8);           % Lab u in cyan
    grid on; xlabel('t [s]'); ylabel('Input [V]'); title('Control Signal','Color','w');
    darken(ax3);

    %f(3)=subplot(3,1,3);
    %plot(t_exp,u)
    %ylabel('Input [V]')
    %xlabel('Time [sec]')
    %grid on

    linkaxes([ax1,ax2,ax3],'x');
end