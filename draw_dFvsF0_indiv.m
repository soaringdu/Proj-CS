function draw_dFvsF0_indiv(time,data_dFvsFo,onset,num_trial,ROI_ind)

    figure('color', 'white', 'position', [300   150   700   500]);   
    plot([onset;onset],repmat(([-1;max(data_dFvsFo(:))]),1,num_trial),'--','color',rgb('grey'),'LineWidth',2);    
    hold on 
    plot(time,data_dFvsFo(:,ROI_ind),'LineWidth',2);    
    xlabel( 'Time  ( s ) ','FontSize',16);
    ylabel( '¦¤F/F0','FontSize',16);
    title(['dFvsF0_raw_trace_of_ROI_',num2str(ROI_ind)],'FontSize',16,'Interpreter','none');
    hold off
    set(gca, 'FontSize',16,'LineWidth',2)
%     saveas(f,[fn_title, '_trace']);    
