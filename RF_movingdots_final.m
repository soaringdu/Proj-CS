
% Functional analysis of rhythmic synaptogenesis project with DXF
    % Original data: 
    %   D:\1_data\节律课题_dxf\original data\WT     # 1-16 
    %   D:\1_data\节律课题_dxf\original data\hcrtr2 # 17-32
    %   D:\1_data\节律课题_dxf\original data\clocka # 33-47
    %   Except 14 15 16; 30,31,32; 45,46,47, use 'end-28:end-21'
    %   Others use 'end-41:end-34'
    %   Batch program. 
    %   Version v20181129 plus F0 calculation by HH method
    %   v20181203 change algorithm for RF size
    %   v20181206 finish codes for statistic analysis
    %   v20190119 add more ROIs by auto-segmentation written by Chen Yang
    %   D:\1_data\节律课题_dxf\original data\RF_movingdots_main_plusHH.m

    clear;clc;close all
    
    %% parameter setting
    rest_dur = 60; % unit:s
    interval = 15; % unit:s; lasting time of a trial from onset to onset
    pre_len = 5; % unit:s; pre stimulus
    post_off_len = 10; % unit:s; post stimulus
    Row_num = 8;
    Col_num = 8;
    iftxt = 1; % if parameter file exists; Y=1 N=0
    Time = 540; % when no paramters file
    ifxls = 0; % if auto-seg needed; Y=1 N=0(import existed results by auto/manual-seg)
    GCaMP_type = 1; % auto-seg; cyto-gcamp=1 nucleus-gcamp=0
    dot_size = 12;
    
    adjrsquare_array = zeros(700,16); 
    RF_size_array = zeros(700,16);
    aver_peak_array = zeros(700,16);
    std_peak_array = zeros(700,16);
    % the number of each group < 16
    % the number of all neurons in one larva < 700
    
for ids = 1:44
        %% import sitmulus seq file
        filename = [num2str(ids),'.mat'];
        rand_seq = load(filename,'index');
        sti_seq = rand_seq.index + 1;

        %% import imaging parameters file
        switch iftxt
            case 1
                data_info = importdata([num2str(ids),'.txt']);
                FrameRate_info = data_info{13};
                if strcmpi(FrameRate_info(end-21:end-14),'interval') == 1
                   Time = str2double(FrameRate_info(end-41:end-34))/1000; 
                   % total time; unit:s
                else
                   Time = str2double(FrameRate_info(end-28:end-21))/1000; 
                   % total time; unit:s
                end
            case 0
                return;
        end

       %% import original calcium data
        switch ifxls
            case 1
                switch GCaMP_type
                    case 1
                        ROI_traces = test_ringfilter2;
                    case 0
                        ROI_traces = test_ringfilter2_n;
                end
                data = ROI_traces'; % col=ROI index
            case 0 
                rawdata = xlsread([num2str(ids),'.xlsx']);
                data = rawdata(:,1:end);
        end
        [num_frame,num_ROI] = size(data);  
        frame_time = Time/num_frame; % unit:s/frame
        time = frame_time*(1:num_frame); % unit:s
            
       %% pre-processing including denoising, detrend, F0 calculation
        dFvsFo = zeros(num_frame,num_ROI);
        smoothwin = 100;
        num_of_std = 3;
        for rr = 1:num_ROI
            data_this_ROI = smooth(data(:,rr),3);
            smoothfactor = round(smoothwin / frame_time);
            baseline_ind_iteration = 1:length(data_this_ROI); 
            baseline_discard = 0;
            count = 0;     
            while ~isempty(baseline_discard) && count <= 100
                baseline_smoothed = smooth(data_this_ROI(baseline_ind_iteration),smoothfactor);
                baseline_smoothed_std = std(baseline_smoothed - data_this_ROI(baseline_ind_iteration));
                baseline_threshold = baseline_smoothed + baseline_smoothed_std * num_of_std;
                baseline_discard = find(data_this_ROI(baseline_ind_iteration) > baseline_threshold);
                baseline_ind_iteration(baseline_discard) = []; 
                peak_ind = setdiff(1:length(data_this_ROI),baseline_ind_iteration); 
                count = count + 1;
            end
            baseline_index = unique([1 baseline_ind_iteration length(data_this_ROI)]);
            data_smooth_baseline = interp1(baseline_index,smooth(baseline_index,data_this_ROI(baseline_index),smoothfactor),1:length(data_this_ROI));
            data_baseline_removed = (data_this_ROI - data_smooth_baseline')./data_smooth_baseline';
            dFvsFo(:,rr) = data_baseline_removed;
            disp(['larvae # ',num2str(ids),'; ROI # ',num2str(rr)]);
        end   
        
        %% segment data
        if ids ~= 36
           num_trial = length(sti_seq) - 1; 
           % discard the last trial since it is not reliable
        else
           num_trial = length(sti_seq) - 2; 
           % specific treatment for clocka #36
        end
        num_trial_real = length(sti_seq); % real sti times
        onset = zeros(1,num_trial); % onset time point of the n trial; unit:s
        pre_sti_frame = zeros(1,num_trial); % unit:frame
        onset_frame = zeros(1,num_trial); % unit:frame
        post_sti_frame = zeros(1,num_trial); % unit:frame
        data_seg = cell(1,num_trial); % segment data by trials

        for ind =1:num_trial
            onset(ind) = rest_dur+(ind-1)*interval;
            pre_sti_frame(ind) = ceil((onset(ind)-pre_len)/frame_time);
            onset_frame(ind) = ceil(onset(ind)/frame_time);
            post_sti_frame(ind) = ceil((onset(ind)+post_off_len)/frame_time);
            data_seg{ind} = dFvsFo(pre_sti_frame(ind):post_sti_frame(ind) ,:);
        end
        
        seg_len = min(post_sti_frame-pre_sti_frame)+1;
        dFvsF0_each_trial = zeros(seg_len,num_ROI,num_trial); 
        % dFvsF0 of each trial of each ROI
        for ind = 1:num_trial
            dFvsF0_each_trial(:,:,ind) = data_seg{ind}(1:seg_len,:);
        end

        %% plot raw traces/heatmap of all/individual ROIs
        draw_dFvsF0_all(time,data_dFvsFo_copy,onset,num_trial);
        draw_dFvsF0_heatmap(time,data_dFvsFo_copy,onset,num_trial,num_ROI,frame_time);
        ROI_ind = 297;
        draw_dFvsF0_indiv(time,dFvsFo,onset,num_trial,ROI_ind);

        %% data storage by directions
        direction_num = Row_num + Col_num;
        index = cell(direction_num,1); 
        dFvsF0_direction = zeros(seg_len,num_ROI,direction_num); 
        last_trial = sti_seq(num_trial+1:end);
        pre_frame = min(onset_frame - pre_sti_frame);
        on_frame = min(onset_frame - pre_sti_frame) + 1; 

        if ids ~= 36
            for ind = 1:direction_num
                index{ind} = find(sti_seq == ind);
                if ind == last_trial(1) 
                   dFvsF0_direction(:,:,ind) = dFvsF0_each_trial(:,:,index{ind}(1));
                else
                   dFvsF0_direction(:,:,ind) = mean(dFvsF0_each_trial(:,:,index{ind}),3); 
                end
            end
        else
            for ind = 1:direction_num
                index{ind} = find(sti_seq == ind);
                if ind == last_trial(1) || ind == last_trial(2) 
                   dFvsF0_direction(:,:,ind) = dFvsF0_each_trial(:,:,index{ind}(1));
                else
                   dFvsF0_direction(:,:,ind) = mean(dFvsF0_each_trial(:,:,index{ind}),3); 
                end
            end
        end

        %% calculate rf
        row_dFvsF0_direction = dFvsF0_direction(:,:,1:8);
        col_dFvsF0_direction = dFvsF0_direction(:,:,9:16);
        win_left = on_frame;
        win_right = round((pre_len+7)/frame_time); 
        [row_peak,~] = max(row_dFvsF0_direction(win_left:win_right,:,:)); % unit:frame
        [col_peak,~] = max(col_dFvsF0_direction(win_left:win_right,:,:)); % unit:frame
        
        res_grid = zeros(Row_num,Col_num,num_ROI);
        for ind = 1:num_ROI
            for i = 1:Row_num
                for j = 1:Col_num
                    res_grid(i,j,ind) = row_peak(1,ind,i) * col_peak(1,ind,j);
                end
            end
        end
        imagesc(res_grid(:,:,297));
        axis square
% % % % % %         sig_row_peak = zeros(num_ROI,8); % evaluate significance
% % % % % %         sig_col_peak = zeros(num_ROI,8); % evaluate significance
% % % % % %         for nn = 1:num_ROI
% % % % % %             for tn = 1:8
% % % % % %                 if row_peak(1,nn,tn) >= mean(row_dFvsF0_direction(1:on_frame-1,nn,tn))+3*std(row_dFvsF0_direction(1:on_frame-1,nn,tn))
% % % % % %                     sig_row_peak(nn,tn) = row_peak(1,nn,tn);
% % % % % %                 else
% % % % % %                     sig_row_peak(nn,tn) = NaN;
% % % % % %                 end
% % % % % %                 if col_peak(1,nn,tn) >= mean(col_dFvsF0_direction(1:on_frame-1,nn,tn))+3*std(col_dFvsF0_direction(1:on_frame-1,nn,tn))
% % % % % %                     sig_col_peak(nn,tn) = col_peak(1,nn,tn);
% % % % % %                 else
% % % % % %                     sig_col_peak(nn,tn) = NaN;
% % % % % %                 end
% % % % % %             end
% % % % % %         end   
% % % % % %  
% % % % % %         tmp = [sig_row_peak sig_col_peak];
% % % % % %         aver_peak_array(1:num_ROI,ids) = mean(tmp,2,'omitnan');
% % % % % %         std_peak_array(1:num_ROI,ids) = std(tmp,0,2,'omitnan');
        
% % % % % %         tmp = zeros(16,1);
% % % % % %         for i = 1:num_ROI
% % % % % %             tmp(1:8,1) = row_peak(1,i,:);
% % % % % %             tmp(9:16,1) = col_peak(1,i,:);
% % % % % %             aver_peak_array(i,ids) = mean(tmp);
% % % % % %             std_peak_array(i,ids) = std(tmp);
% % % % % %         end
        
       %% 2D Gaussian Fitting
        x_range = 1:8;
        y_range = 1:8;
        [Y,X] = meshgrid(x_range,y_range);
        % c1*exp(-((x-a).^2/d1^2+(y-b).^2/d2^2)) 
        % c0 + c1 * exp(-((x-c2)*cos(c3)+(y-c4)*sin(c3)).^2/(2*c5^2) - ( - (x-c2)*sin(c3)+(y-c4)*cos(c3)).^2 /(2*c6^2)); 
        % c2和c4的初值很重要，[1,8]，复杂版本含旋转的二维高斯拟合函数
        % c2和c4对应中心坐标，重要
        % c5和c6对应两个方向上的sigma，重要
        % c3对应旋转角度
        coeff = cell(1,num_ROI); % 储存全部的拟合参数
        gof = cell(1,num_ROI); % 储存全部的拟合好坏判定结果
        gof_each = zeros(1,num_ROI); % r-square，越接近1，拟合质量越好
        for i = 1:num_ROI
            [coeff{i},gof{i}] = fit_gauss2d(X, Y, res_grid(:,:,i));
            gof_each(1,i) = gof{i}.adjrsquare;
            disp(['larvae # ',num2str(ids),'; ROI # ',num2str(i)]);
        end
end
% xlswrite('D:\1_data\节律课题_dxf\original data\20190119_final\aver_sig_peak_hcrtr2_190330.xlsx',aver_peak_array); 
% xlswrite('D:\1_data\节律课题_dxf\original data\20190119_final\std_sig_peak_hcrtr2_190330.xlsx',std_peak_array);

    %% plot averaged trace of one direction of one ROI
    ROI_ind = 25;
    Direction_ind = 8;
    tmp = dFvsF0_direction(:,ROI_ind,Direction_ind); 
    plot(tmp);
    hold on 
    x1 = [on_frame on_frame];
    y1 = [min(tmp) max(tmp)];
    line(x1,y1,'LineWidth',3);
    xlabel( 'frames ','FontSize',16);
    ylabel( 'ΔF/F0','FontSize',16);
    title(['ROI_' num2str(ROI_ind) '_Direction_' num2str(Direction_ind)],'FontSize',16,'Interpreter','none');
    box off
    hold off
    set(gca, 'FontSize',16,'LineWidth',2)
    
    %% plot traces of one direction of one ROI
    ROI_ind = 101;
    Direction_ind = 5;
    tmp = F_direction_each_trial(:,ROI_ind,Direction_ind); 
    plot(tmp);
    hold on 
    x1 = [on_frame on_frame];
    y1 = [min(tmp) max(tmp)];
    line(x1,y1,'LineWidth',3);
    x2 = [on_frame+seg_len on_frame+seg_len];
    y2 = [min(tmp) max(tmp)];
    line(x2,y2,'LineWidth',3);
    xlabel( 'frames ','FontSize',16);
    ylabel( 'ΔF/F0','FontSize',16);
    title(['ROI_' num2str(ROI_ind) '_Direction_' num2str(Direction_ind)],'FontSize',16,'Interpreter','none');
    box off
    hold off
    set(gca, 'FontSize',16,'LineWidth',2)
    
    %% plot responses to all directions of one ROI
    ROI_ind = 297;   
    tmp1 = zeros(direction_num,seg_len); % mean response
    tmp2 = zeros(direction_num,seg_len*2); % each trial
    for ind = 1:direction_num
        tmp1(ind,:) = dFvsF0_direction(:,ROI_ind,ind);
        tmp2(ind,:) = F_direction_each_trial(:,ROI_ind,ind);
    end
    imagesc(tmp2);    
    x1 = [on_frame on_frame];
    y1 = [1 direction_num];
    line(x1,y1,'LineWidth',3,'color',rgb('white')); % line at the first frame after sti
    x2 = [on_frame+seg_len on_frame+seg_len];
    y2 = [1 direction_num];
    line(x2,y2,'LineWidth',3,'color',rgb('white')); % line at the first frame after sti
    xlabel_frame = 0:2:seg_len*2;
    ylabel_direction = 1:1:direction_num;
    xlabel( 'frames','FontSize',16);
    ylabel( 'Direction','FontSize',16);
    set(gca,'XTick',xlabel_frame,'YTick',ylabel_direction, 'FontSize',16,'LineWidth',2)
    
    %% plot responses of all ROIs to one direction
    Direction_ind = 6;   
    tmp2 = zeros(num_ROI,seg_len);
    for ind = 1:num_ROI
        tmp2(ind,:) = dFvsF0_direction(:,ind,Direction_ind);
    end
    imagesc(tmp2);
    x2 = [on_frame on_frame];
    y2 = [1 num_ROI];
    line(x2,y2,'LineWidth',3,'color',rgb('white'));
    xlabel_frame = 1:1:seg_len;
    ylabel_roi = 1:1:num_ROI;
    xlabel( 'frames','FontSize',16);
    ylabel( 'ROI','FontSize',16);
    set(gca,'XTick',xlabel_frame,'YTick',ylabel_roi, 'FontSize',16,'LineWidth',2)
    
   %% plot averaged responses of all ROIs to all directions
    ROI_ind = 297;
    x1 = [on_frame on_frame];
    y1 = [0 direction_num+1];
    tmp3 = zeros(seg_len,direction_num);
    for i = 1:direction_num   
        tmp3(:,i) = dFvsF0_direction(:,ROI_ind,i);
    end
    imagesc(tmp3');
    line(x1,y1,'LineWidth',3,'color',rgb('white'),'LineStyle','--');
    xlabel_time = 0:3:15;
    xlabel_frame = round(xlabel_time/frame_time);
    set(gca,'XTick',xlabel_frame,'XTicklabel',xlabel_time, 'FontSize',16,'LineWidth',2);
    box off
    xlabel( 'Time (s)','FontSize',16);
    ylabel( 'Direction','FontSize',16);

    %% investigate the distribution of r-square,20180722
    fn_1 = uigetfile('*.xlsx','Select the related xlsx file'); 
    cumulative_result = xlsread(fn_1);
    edge = 0:0.1:1;
    wt_cum = cumulative_result(1:16,:);
    hcr_cum = cumulative_result(17:32,:);
    clock_cum = cumulative_result(33:47,:);
    figure
    plot(edge,mean(wt_cum),'k','LineWidth',4);
    hold on 
    plot(edge,mean(hcr_cum),'r','LineWidth',4);
    plot(edge,mean(clock_cum),'b','LineWidth',4);
    hold off
    legend('wt','hcrtr2','clocka');
    ylabel('cumulative percentage');
    xlabel('r-square');
    
    figure
    p1 = shadedErrorBar(edge,mean(wt_cum),std(wt_cum),'lineprops','-k','patchSaturation',0.2,'transparent',true);
    hold on 
    p2 = shadedErrorBar(edge,mean(hcr_cum),std(hcr_cum),'lineprops','-r','patchSaturation',0.2,'transparent',true);
    p3 = shadedErrorBar(edge,mean(clock_cum),std(clock_cum),'lineprops','-b','patchSaturation',0.2,'transparent',true);
    line([0.6 0.6],[0 1],'Color','k','LineStyle','--');
    hold off
    legend('wt','','hcrtr2','','clocka');
    ylabel('cumulative percentage');
    xlabel('r-square');
    ylim([0 1]);
    xlim([0 1]);
    grid off
    box off
   
    a = mean(wt_cum);
    b = mean(hcr_cum);
    c = mean(clock_cum);
    [h,p] = kstest2(sib,mut); % 对sib、mut两组数据做假设检验，h = 1，则两组数没有差别。h = 0，则两组数有差异。
    if h == 1
        disp(['两组数据有显著性差异，p =',num2str(p)]);
    else
        disp('两组数据没有差异');
    end

    %% compare RF of different groups with different thresholds of adjrsquare
    num_ROI_each = [1,20,33;12,29,44;12,10,12];
    filename_wt_rf = uigetfile('*.xlsx','Select the data file'); % select the excel file
    filename_hcr_rf = uigetfile('*.xlsx','Select the data file'); % select the excel file
    filename_clock_rf = uigetfile('*.xlsx','Select the data file'); % select the excel file
    rf = cell(1,3);
    [rf{1}, ~, ~] = xlsread(filename_wt_rf);
    [rf{2}, ~, ~] = xlsread(filename_hcr_rf);
    [rf{3}, ~, ~] = xlsread(filename_clock_rf);
    rf_all = cell(1,3);
    rf_all{1} = rf{1}(rf{1}~=0); % 1-wt,2-hcr,3-clock
    rf_all{2} = rf{2}(rf{2}~=0);
    rf_all{3} = rf{3}(rf{3}~=0);
    
    rf_each_group = cell(1,3);
    for ind_group = 1:3 
        tmp_group = rf{ind_group}(:,num_ROI_each(1,ind_group):num_ROI_each(2,ind_group));
        rf_each_group{ind_group} = cell(1,num_ROI_each(3,ind_group));
        for ind_larva = 1:num_ROI_each(3,ind_group)
            tmp_larva = tmp_group(:,ind_larva);
            tmp_larva = tmp_larva(tmp_larva~=0); % rf of all neurons in current larva
            rf_each_group{ind_group}{ind_larva} = tmp_larva;
        end
    end
    
    filename_wt_adjr = uigetfile('*.xlsx','Select the data file'); % select the excel file
    filename_hcr_adjr = uigetfile('*.xlsx','Select the data file'); % select the excel file
    filename_clock_adjr = uigetfile('*.xlsx','Select the data file'); % select the excel file
    adjr = cell(1,3);
    [adjr{1}, ~, ~] = xlsread(filename_wt_adjr);
    [adjr{2}, ~, ~] = xlsread(filename_hcr_adjr);
    [adjr{3}, ~, ~] = xlsread(filename_clock_adjr);
    adjr_all = cell(1,3);
    adjr_all{1} = adjr{1}(adjr{1}~=0); 
    adjr_all{2} = adjr{2}(adjr{2}~=0);
    adjr_all{3} = adjr{3}(adjr{3}~=0);
    
    adjr_each_group = cell(1,3);
    for ind_group = 1:3 
        tmp_group = adjr{ind_group}(:,num_ROI_each(1,ind_group):num_ROI_each(2,ind_group));
        adjr_each_group{ind_group} = cell(1,num_ROI_each(3,ind_group));
        for ind_larva = 1:num_ROI_each(3,ind_group)
            tmp_larva = tmp_group(:,ind_larva);
            tmp_larva = tmp_larva(tmp_larva~=0); % rf of all neurons in current larva
            adjr_each_group{ind_group}{ind_larva} = tmp_larva;
        end
    end
    
    %% all adjr > 0, cumulative under diff adjrsquare
    bin = 25;
    xbin = 6:6:150;
    aver_rf = zeros(3,10);
    aver_tmp_rf = zeros(10,3);
    std_tmp_rf = zeros(10,3);
    sem_tmp_rf = zeros(10,3);
    location = zeros(3,bin);
    fraction = zeros(3,bin);
    for i = 1:10
        threshold = (i-1)*0.1;
        subplot(2,5,i);
        roi_n = zeros(1,3);
        tmp_ks = cell(1,3); % for kstest2
        for ind_group = 1:3 % 1-wt,2-hcr,3-clock    
            ind_rf = find(adjr_all{ind_group} > threshold);
            tmp_rf = rf_all{ind_group}(ind_rf);
            aver_tmp_rf(i,ind_group) = mean(tmp_rf);
            std_tmp_rf(i,ind_group) = std(tmp_rf);
            tmp_ks{ind_group} = tmp_rf;
            aver_rf(ind_group,i) = mean(tmp_rf);
            [n,b] = hist(tmp_rf,xbin);
            location(ind_group,:) = b;
            fraction(ind_group,:) = cumsum(n)/sum(n);
            plot(b,cumsum(n)/sum(n),'LineWidth',3);
            hold on
            roi_n(ind_group) = length(tmp_rf);
            sem_tmp_rf(i,ind_group) = std_tmp_rf(i,ind_group)/sqrt(roi_n(ind_group));
        end
        [h_hw,p_hw] = kstest2(tmp_ks{2},tmp_ks{1},'Tail','Smaller');
        [h_cw,p_cw] = kstest2(tmp_ks{3},tmp_ks{1},'Tail','Smaller');
        [h_ch,p_ch] = kstest2(tmp_ks{3},tmp_ks{2},'Tail','Smaller');
        legend(['wt n = ',num2str(roi_n(1))],['hcrtr2 n = ',num2str(roi_n(2))],['clocka n = ',num2str(roi_n(3))],'location','southeast');
        legend('boxoff')
        title({['adjsquare > ',num2str(threshold)];['h-hw=',num2str(h_hw),' h-cw=',num2str(h_cw),' h-ch=',num2str(h_ch)];['p-hw=',num2str(p_hw)];['p-cw=',num2str(p_cw)];['p-ch=',num2str(p_ch)]});
        xlabel('RF size /degree');
        ylabel('Percentage / %');
        box off
    end
   
    
    
    
    
    
    

