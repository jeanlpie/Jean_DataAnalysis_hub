%% Plot several firing frequencies in one plot

fc = load('z-score_PPC_28_T_R9_Co_Off')
red_line = fc.Zpsth

fc = load('z-score_PPC_28_V_R4_Co_Off')
blue_line = fc.Zpsth

fc = load('z-score_PPC_28_M_R13_Co_Off.mat')
purple_line = fc.Zpsth

labels = {'TR', 'VR', 'MR'} % Red blue purple labels
filename=strcat('frequencies_newplot.pdf');

%Make SDF with Z-scores
figure()
[B] = msdf([((edges+BinSize*1e3/2)*1e-6)', red_line'],'Gauss',2);
plot(B(:,1),B(:,2), 'r')
hold on
[B] = msdf([((edges+BinSize*1e3/2)*1e-6)', blue_line'],'Gauss',2);
plot(B(:,1),B(:,2), 'b')
[B] = msdf([((edges+BinSize*1e3/2)*1e-6)', purple_line'],'Gauss',2);
plot(B(:,1),B(:,2), 'color', [0.49 0.18 0.45])

% title('Correct')
legend(labels, 'Location', 'NorthWest')
xlabel('Time (sec)') %Label x-axis
ylabel('psth z-score') %Label y-axis
set(findall(gcf,'-property','FontSize'),'FontSize',18)
xlim([win(1) win(2)]*1e-3);
ylim([-1 1]*1.5);

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf,filename,'-dpdf')
pause(0.1)



% boolBslnNorm    = 0;
% boolZscored     = 0;
% boolGaussFilt   = 1;
% FilterSize      = 25; %in ms