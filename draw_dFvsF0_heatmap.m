function draw_dFvsF0_heatmap(time,data_dFvsFo,onset,num_trial,num_ROI,frame_time)

%     f_heatmap = figure;
    clims = [0 1];
    imagesc(data_dFvsFo',clims);
    hold on
    plot([onset/frame_time;onset/frame_time],repmat(([1;num_ROI]),1,num_trial),'--','color',rgb('white'),'LineWidth',1);     
    xlabel( 'Time  ( s ) ','FontSize',16);
    ylabel( 'Cell No.','FontSize',16);
    hold off   
    xlabel_time = 0:60:time(end);
    xlabel_frame = round(xlabel_time/frame_time);
    set(gca,'XTick',xlabel_frame,'XTicklabel',xlabel_time, 'FontSize',16,'LineWidth',2)
    % saveas(f_heatmap,[fn_title, '_heatmap']);