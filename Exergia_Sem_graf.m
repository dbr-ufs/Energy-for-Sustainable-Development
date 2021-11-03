function Exergia_Sem_graf(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 06-Jan-2018 16:27:07

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.12598393574297 0.415699481865286 0.462369477911645 0.542849740932642]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'MarkerSize',8,'Marker','square','Color',[0 0.5 0]);
set(plot1(2),'MarkerSize',8,'Marker','square','Color',[0 0 1]);
%set(plot1(3),'Marker','o','Color',[0 0.5 0]);
set(plot1(3),'Visible','off','Marker','o','Color',[0 0.5 0]);
set(plot1(4),'Marker','o','Color',[0 0 1]);
set(plot1(5),'Marker','*','Color',[0 0.5 0]);
set(plot1(6),'Marker','*','Color',[0 0 1]);
%set(plot1(7),'MarkerSize',8,'Marker','x','Color',[0 0.5 0]);
set(plot1(7),'Visible','off','MarkerSize',8,'Marker','x','Color',[0 0.5 0]);
set(plot1(8),'MarkerSize',8,'Marker','x','Color',[0 0 1]);
set(plot1(9),'Marker','diamond','Color',[0 0.5 0]);
set(plot1(10) ,'Marker','diamond','Color',[0 0 1]);
set(plot1(11),'Color',[0.400000005960464 1 0.400000005960464],'LineWidth',2);
set(plot1(12),'Color',[0 0.600000023841858 1],'LineWidth',2);

% Create xlabel
xlabel('Cell Temperature [�C]','FontName','Times New Roman');

% Create ylabel
ylabel('Exergetic Efficiency [%]','FontName','Times New Roman');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[30 110]);
% Uncomment the following line to preserve the Y-limits of the axes
%ylim(axes1,[3 13]);
ylim(axes1,[4.8 15]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',11,'XGrid','on','YGrid',...
    'on');
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.118331777815488 0.50636422058883 0.157515777395298 0.0777202072538863],...
%     'String','\eta_{ref} = 5%',...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.116816865045694 0.777859449222783 0.164873275799379 0.0785714269393967],...
%     'String','\eta_{ref} = 10%',...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

%%%% ABAIXO CALOR

% Create textbox
annotation(figure1,'textbox',...
    [0.176681878419111 0.418508923431208 0.157515777395298 0.0777202072538863],...
    'String','\eta_{ref} = 5%',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.173154893214708 0.68483619340883 0.164873275799379 0.0785714269393967],...
    'String','\eta_{ref} = 10%',...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');
