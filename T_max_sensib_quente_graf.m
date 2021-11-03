function T_max_sensib_quente_graf(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 07-Jan-2018 17:55:26

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.12598393574297 0.415699481865286 0.462369477911645 0.542849740932642]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Color',[0 0.498039215803146 0],...
    'Parent',axes1);
set(plot1(1),'Marker','+');
set(plot1(2),'Marker','o');
set(plot1(3),'Visible','off','Marker','*');
set(plot1(4),'Marker','^');
set(plot1(5),'Marker','x');
set(plot1(6),'Visible','off','Marker','square');
set(plot1(7),'Visible','off','Marker','diamond');
set(plot1(8),'Visible','off','Marker','v');

% Create xlabel
xlabel('Solar Irradiance [W/m�]','FontName','Times New Roman');

% Create ylabel
ylabel('Sensitivity [%]','FontName','Times New Roman');

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-100 100]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',11,'XGrid','on','YGrid',...
    'on');
% Create textbox
annotation(figure1,'textbox',...
    [0.402724137931033 0.828875417711564 0.140086203363949 0.0913705564089837],...
    'String','F_R.U_L',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.358142857142855 0.457834058933268 0.183189650355228 0.0913705564089838],...
    'String','F_R.(\tau.\alpha)_e',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.530556650246303 0.474654325031285 0.07543103287703 0.0736040594765379],...
    'String','G',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.486775862068964 0.680266843651596 0.120689652217873 0.0913705564089838],...
    'String','\eta_{ref}',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');
