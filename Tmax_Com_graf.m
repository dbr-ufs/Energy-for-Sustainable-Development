function Tmax_Com_graf(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 20-Dec-2018 21:38:20

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.12598393574297 0.415699481865286 0.462369477911645 0.542849740932642]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'Marker','square','Color',[0 0.498039215803146 0]);
set(plot1(2),'Visible','off');
set(plot1(3),'MarkerSize',8,'Marker','x','Color',[0 0 1]);
set(plot1(4),'Visible','off','MarkerSize',8,'Marker','x','Color',[1 0 0]);
set(plot1(5),'Marker','diamond','Color',[0 0 1]);
set(plot1(6),'Marker','*','Color',[0 0.498039215803146 0]);
set(plot1(7),'Marker','*','Color',[0 0 1]);
set(plot1(8),'Marker','*','Color',[1 0 0]);
set(plot1(9),'Visible','off');
set(plot1(10),'Marker','x','Color',[1 0 0]);
set(plot1(11),'Visible','off');
set(plot1(12),'Marker','diamond','Color',[0 0 1]);
set(plot1(13),'Marker','diamond','Color',[1 0 0]);

% Create xlabel
xlabel('Solar Irradiance [W/m^2]','FontName','Times New Roman');

% Create ylabel
ylabel('Maximum Temperature [�C]','FontName','Times New Roman');

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 400]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',11,'XGrid','on','YGrid',...
    'on');
% Create textbox
annotation(figure1,'textbox',...
    [0.130854918503049 0.888095238095239 0.117359367211237 0.0572326277816494],...
    'String','\eta_{ref} = 5%',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create arrow
annotation(figure1,'arrow',[0.204556986729117 0.226785714285714],...
    [0.888920555753474 0.857142857142857],'HeadWidth',6,'HeadLength',6);

% Create textbox
annotation(figure1,'textbox',...
    [0.379377698146212 0.816666666666668 0.127765158996644 0.0589611894271268],...
    'String','\eta_{ref} = 20%',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.286338915380487 0.854761904761906 0.0600896560480836 0.0509922762794542],...
    'String','10%',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

