%% Coefficient and Simulation

N_repeat = 5; % Number of different E repeats
N_repeat2 = 3; % Number of repeats for the same E

% Record all the data during repeat
burst_freq_list_sum = []; % Burst frequency for all repeats
burst_size_list_sum = []; % Volume integral with time for all repeats
burst_waittime_list_sum = []; % Time until burst for all repeats
burst_duringtime_list_sum = []; % Burst duration for all repeats

% Repeat for generating different E distributions
for repeat_i = 1:N_repeat

N = 10000; % Number of nucleation sites
dt_rec = 1; % Time interval of recording data

% Parameter for E* = 0.07 of lognormal distribution
sigma_lgn = 0.4;
mu_lgn = -1.6;
Y = 0;
E = exp(mu_lgn + sigma_lgn.*randn(N,1));
E = sort(E);

% E_max = 1;
% theta = 2; %x^theta distribution
% Y = 0;
% E = E_max*rand(N,1).^(1/(theta+1));
% E = sort(E);

kRn = 5; % R_n = kRn*R_c
noise_ratio = 0.0; % Noise ratio for R_n, 0 means no noise, 1 means [0,2*kRn*R_c]

k_dis = 1e-7; % Dissolution rate parameter

growth_num_max = 500; % Number of maximum growth times

% Repeat for one E distribution
for repeat2_i = 1:N_repeat2

tic;
run nucleation_GrowthTime_pipesimulation.m;
toc;

%% Burst Calculation
growth_num_thresh = 3; % Minimum growth times for calculation
N_growth_max = find(growth_num_sum >= growth_num_thresh,1,'last');

growthsum_count = 0;
burst_count = 0;
burst_idx_list = [];
growth_duringtime_list = [];

burst_check_sum = 0;
for idx_i = N_growth_max:-1:1
    if growth_num_sum(idx_i) < growth_num_thresh
        continue
    end
    growthsum_count = growthsum_count + 1;

    t_growth_list_i = t_growth_list_sum{idx_i};
    t_dissolution_list_i = t_dissolution_list_sum{idx_i};
            
    t_delay_list = [];
    t_during_list = [];
    for ti = 1:length(t_growth_list_i)
        t_growth = t_growth_list_i(ti);
        t_dissolution = t_dissolution_list_i(find(t_dissolution_list_i<=t_growth,1,'last'));
        t_dissolution_after = t_dissolution_list_i(find(t_dissolution_list_i>t_growth,1,'first'));
        if isempty(t_dissolution)
            continue;
        end
        if isempty(t_dissolution_after)
            continue;
        end
        t_delay = t_growth-t_dissolution;
        t_delay_list = [t_delay_list,t_delay];

        t_during = t_dissolution_after - t_growth;
        t_during_list = [t_during_list,t_during];
    end

    if isempty(t_delay_list)
        t_delay_list = 0;
    end
    if isempty(t_during_list)
        t_during_list = 0;
    end
    growth_duringtime_list = [growth_duringtime_list,mean(t_during_list)]; % Record lifetimes for all condnesate

    if mean(t_delay_list)/mean(t_during_list) >= 1 && burst_check_sum < 5 % Burst condensates have longer time until burst than burst duration
        burst_count = burst_count + 1;
        burst_idx_list = [burst_idx_list,idx_i];
    else
        burst_check_sum = burst_check_sum + 1;
    end

end
burst_fraction = burst_count/sum(growth_num_sum>=growth_num_thresh)

%% Burst Frequency and Size Calculation
burst_freq_list = zeros(burst_count,1); % Burst frequency
burst_size_list = zeros(burst_count,1); % Volume integral with time
burst_waittime_list = zeros(burst_count,1); % Time until burst
burst_duringtime_list = zeros(burst_count,1); % Burst duration

t_delay_list_sum = cell(burst_count,1);

for burst_i = 1:burst_count
    burst_idx = burst_idx_list(burst_i);
    
    t_growth_list_i = t_growth_list_sum{burst_idx};
    t_dissolution_list_i = t_dissolution_list_sum{burst_idx};
    V_final_i = V_final_sum{burst_idx};
    if isempty(V_final_i)
        V_final_i = 0;
    end
    V_intgrowth_i = V_intgrowth_sum{burst_idx};
    if isempty(V_intgrowth_i)
        V_intgrowth_i = 0;
    end

    t_delay_list = [];
    t_during_list = [];
    burst_size_i = [];
    for ti = 1:length(t_growth_list_i)
        t_growth = t_growth_list_i(ti);
        t_dissolution = t_dissolution_list_i(find(t_dissolution_list_i<=t_growth,1,'last'));
        t_dissolution_after = t_dissolution_list_i(find(t_dissolution_list_i>t_growth,1,'first'));
        if isempty(t_dissolution)
            continue;
        end
        if isempty(t_dissolution_after)
            continue;
        end
        t_delay = t_growth-t_dissolution;
        t_delay_list = [t_delay_list,t_delay];

        t_during = t_dissolution_after - t_growth;
        t_during_list = [t_during_list,t_during];

        % t_life = t_dissolution_after - t_growth;
        % burst_size = t_life * V_final_i(t_dissolution_list_i==t_dissolution_after) / 2;
        % burst_size_i = [burst_size_i,burst_size];

        burst_size_i = [burst_size_i, V_intgrowth_i(t_dissolution_list_i==t_dissolution_after)];

    end
    if isempty(t_delay_list)
        t_delay_list = 0;
    end
    if isempty(t_during_list)
        t_during_list = 0;
    end
    if isempty(burst_size_i)
        burst_size_i = 0;
    end

    t_delay_list_sum{burst_i} = t_delay_list;

    t_growth_diff_list_i = diff(t_growth_list_i);

    burst_freq_list(burst_i) = 1./mean(t_delay_list);
    burst_size_list(burst_i) = mean(burst_size_i);
    burst_waittime_list(burst_i) = mean(t_delay_list);
    burst_duringtime_list(burst_i) = mean(t_during_list);

end


%% Simulation Data Distribution Calculation for Repeat

burst_freq_list_sum = [burst_freq_list_sum;burst_freq_list];
burst_size_list_sum = [burst_size_list_sum;burst_size_list];
burst_waittime_list_sum = [burst_waittime_list_sum;burst_waittime_list];
burst_duringtime_list_sum = [burst_duringtime_list_sum;burst_duringtime_list];

end

end



%% Read Experimental Data From Table

burst_data_table = readtable('./CAST_burst_data.csv');

[N_gene,~] = size(burst_data_table);

k_on_list = zeros(N_gene,1);
k_off_list = zeros(N_gene,1);
burst_data_freq_list = zeros(N_gene,1);
burst_data_size_list = zeros(N_gene,1);

for gene_idx = 1:N_gene
    temp_list = str2num(cell2mat(burst_data_table{gene_idx,2}));
    k_on_list(gene_idx) = temp_list(1);
    k_off_list(gene_idx) = temp_list(2);

    temp_list = str2num(cell2mat(burst_data_table{gene_idx,3}));
    burst_data_freq_list(gene_idx) = temp_list(1);

    temp_list = str2num(cell2mat(burst_data_table{gene_idx,4}));
    burst_data_size_list(gene_idx) = temp_list(1);
end


%% Plot Patch Figures for Comparison

t_decay = 3.9/log(2); % 1/mRNA degradation rate in the unit hour



% log10 burst frequency
figure;
color_list = [0.01,0.29,0.63;0.76,0.22,0.19];

[counts,edges] = histcounts(log10(k_on_list/t_decay),-3.5:0.15:0,'Normalization','PDF');
edges_plot = edges(1:end-1) + (edges(2)-edges(1))/2;
patch([edges_plot(1),edges_plot,edges_plot(end)],[0,counts,0],color_list(1,:),'FaceAlpha',0.1,'EdgeAlpha',0);
hold on
f1 = plot(edges_plot,counts,'Color',color_list(1,:),'LineWidth',1.5);
counts_exp = counts*0.15;
counts_exp(counts_exp==0) = 1e-3;

[counts,edges] = histcounts(log10(1./burst_waittime_list_sum/(mean(1./burst_waittime_list_sum))*mean(k_on_list/t_decay)),-3.5:0.15:0,'Normalization','PDF');
edges_plot = edges(1:end-1) + (edges(2)-edges(1))/2;
patch([edges_plot(1),edges_plot,edges_plot(end)],[0,counts,0],color_list(2,:),'FaceAlpha',0.1,'EdgeAlpha',0);
hold on
f2 = plot(edges_plot,counts,'Color',color_list(2,:),'LineWidth',1.5);
counts_simul = counts*0.15;
counts_simul(counts_simul==0) = 1e-3;

set(gca,'FontSize',20)
set(gca,'Box','on')
xlabel('log_{10}(burst frequency) (1/h)','Fontsize',25)
ylabel('Probability Density','Fontsize',25)
legend([f1,f2],'Experiment','Simulation')

disp(sum(counts_exp.*log(counts_exp./counts_simul))); % K-L Divergence for burst frequency





% Time until burst
figure;
color_list = [0.01,0.29,0.63;0.76,0.22,0.19];

[counts,edges] = histcounts(t_decay./k_on_list,0:2:100,'Normalization','PDF');
edges_plot = edges(1:end-1) + (edges(2)-edges(1))/2;
patch([edges_plot(1),edges_plot,edges_plot(end)],[0,counts,0],color_list(1,:),'FaceAlpha',0.1,'EdgeAlpha',0);
hold on
f1 = plot(edges_plot,counts,'Color',color_list(1,:),'LineWidth',1.5);
counts_exp = counts*2;
counts_exp(counts_exp==0) = 1e-3;

[counts,edges] = histcounts(burst_waittime_list_sum.*(mean(1./burst_waittime_list_sum))/mean(k_on_list/t_decay),0:2:100,'Normalization','PDF');
edges_plot = edges(1:end-1) + (edges(2)-edges(1))/2;
patch([edges_plot(1),edges_plot,edges_plot(end)],[0,counts,0],color_list(2,:),'FaceAlpha',0.1,'EdgeAlpha',0);
hold on
f2 = plot(edges_plot,counts,'Color',color_list(2,:),'LineWidth',1.5);
counts_simul = counts*2;
counts_simul(counts_simul==0) = 1e-3;

set(gca,'FontSize',20)
set(gca,'Box','on')
xlabel('Time until burst (h)','Fontsize',25)
ylabel('Probability Density','Fontsize',25)
legend([f1,f2],'Experiment','Simulation')

disp(sum(counts_exp.*log(counts_exp./counts_simul))); % K-L Divergence for time until burst





% Burst duration
figure;
color_list = [0.01,0.29,0.63;0.76,0.22,0.19];

[counts,edges] = histcounts(t_decay./k_off_list,0:0.2:10,'Normalization','PDF');
edges_plot = edges(1:end-1) + (edges(2)-edges(1))/2;
patch([edges_plot(1),edges_plot,edges_plot(end)],[0,counts,0],color_list(1,:),'FaceAlpha',0.1,'EdgeAlpha',0);
hold on
f1 = plot(edges_plot,counts,'Color',color_list(1,:),'LineWidth',1.5);
counts_exp = counts*0.2;
counts_exp(counts_exp==0) = 1e-3;

[counts,edges] = histcounts(burst_duringtime_list_sum.*(mean(1./burst_waittime_list_sum))/mean(k_on_list/t_decay),0:0.2:10,'Normalization','PDF');
edges_plot = edges(1:end-1) + (edges(2)-edges(1))/2;
patch([edges_plot(1),edges_plot,edges_plot(end)],[0,counts,0],color_list(2,:),'FaceAlpha',0.1,'EdgeAlpha',0);
hold on
f2 = plot(edges_plot,counts,'Color',color_list(2,:),'LineWidth',1.5);
counts_simul = counts*0.2;
counts_simul(counts_simul==0) = 1e-3;

set(gca,'FontSize',20)
set(gca,'Box','on')
xlabel('Burst duration (h)','Fontsize',25)
ylabel('Probability Density','Fontsize',25)
legend([f1,f2],'Experiment','Simulation')
disp(sum(counts_exp.*log(counts_exp./counts_simul))); % K-L Divergence for burst duration



