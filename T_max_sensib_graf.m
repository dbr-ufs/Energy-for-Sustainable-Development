function T_max_sensib_graf(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 07-Jan-2018 16:45:21

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.12598393574297 0.415699481865286 0.462369477911645 0.542849740932642]);
ylim(axes1,[-40 110]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Color',[0 0.498039215803146 0],...
    'Parent',axes1);
set(plot1(1),'Visible','off','Marker','+');
set(plot1(2),'Marker','o');
set(plot1(3),'Visible','off','Marker','*');
set(plot1(4),'Marker','^');
set(plot1(5),'Marker','x');
set(plot1(6),'Marker','square');
set(plot1(7),'Marker','diamond');
set(plot1(8),'Marker','v');

% Create xlabel
xlabel('Solar Irradiance [W/m�]');

% Create ylabel
ylabel('Sensitivity [%]');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',11,'XGrid','on','YGrid',...
    'on');
% Create textbox
annotation(figure1,'textbox',...
    [0.504017241379309 0.862545324664254 0.120689652217873 0.0913705564089838],...
    'String','T_{amb}',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.463746305418717 0.673023932350989 0.183189650355228 0.0913705564089837],...
    'String','F_R.(\tau.\alpha)_e',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.517625615763544 0.6265530592085 0.0754310328770301 0.0736040594765378],...
    'String','G',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.474214285714285 0.506006770156158 0.155172409810897 0.0913705564089837],...
    'String','T_{cooling}',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.517440886699506 0.40848441042919 0.118534479868309 0.0736040594765379],...
    'String','COP',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.232465517241378 0.401027316445741 0.140086203363949 0.0913705564089838],...
    'String','F_R.U_L',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

