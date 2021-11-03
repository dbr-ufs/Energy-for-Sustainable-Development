function Tmax_Sem_frio_graf(X1, YMatrix1)
%CREATEFIGURE1(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 09-Jan-2018 18:14:27

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.12598393574297 0.415699481865286 0.462369477911645 0.542849740932642]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1,'Parent',axes1);
set(plot1(1),'Marker','square','Color',[0 0.498039215803146 0]);
set(plot1(2),'Visible','off','Marker','square','Color',[0 0 1]);
set(plot1(3),'Visible','off','Marker','square','Color',[1 0 0]);
set(plot1(4),'Marker','o','Color',[0 0.498039215803146 0]);
set(plot1(5),'Marker','o','Color',[0 0 1]);
set(plot1(6),'Visible','off','Marker','o','Color',[1 0 0]);
set(plot1(7),'Marker','*','Color',[0 0.498039215803146 0]);
set(plot1(8),'Marker','*','Color',[0 0 1]);
set(plot1(9),'Marker','*','Color',[1 0 0]);
set(plot1(10),'Marker','x','Color',[0 0.498039215803146 0]);
set(plot1(11),'Marker','x','Color',[0 0 1]);
set(plot1(12),'Marker','x','Color',[1 0 0]);
set(plot1(13),'Visible','off','Marker','diamond',...
    'Color',[0 0.498039215803146 0]);
set(plot1(14),'Marker','diamond','Color',[0 0 1]);
set(plot1(15),'Marker','diamond','Color',[1 0 0]);

% Create xlabel
xlabel('Solar Irradiance [W/m^2]','FontName','Times New Roman');

% Create ylabel
ylabel('Maximum Temperature [�C]','FontName','Times New Roman');

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[40 120]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',11,'XGrid','on','YGrid',...
    'on');
% Create textbox
annotation(figure1,'textbox',...
    [0.425800492610837 0.622166246851385 0.153940886699508 0.0478589420654912],...
    'String',{'\eta_{ref} = 10%'},...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.424568965517241 0.430730478589421 0.157327586206896 0.0453220582943507],...
    'String',{'\eta_{ref} = 20%'},...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textarrow
annotation(figure1,'textarrow',[0.379310344827587 0.413793103448276],...
    [0.860460957178841 0.823677581863979],'TextBackgroundColor',[1 1 1],...
    'String',{'\eta_{ref} = 5%'},...
    'HeadWidth',6,...
    'HeadLength',6,...
    'FontName','Times New Roman');

