%% 4.1
load('C:\Users\Soren\Desktop\Dataset\AD.mat');  
load('C:\Users\Soren\Desktop\Dataset\Normal.mat');  

% Set parameters
sampling_rate = 200;
frequency_range = [35 40];

%set initial variables
% calculate number of AD patients
number_AD = numel(AD);
% calculate number of Normal patients
num_Normal = numel(normal); 
plv_rare_AD = zeros(1, number_AD);
plv_frequent_AD = zeros(1, number_AD);
plv_rare_Normal = zeros(1, num_Normal);
plv_frequent_Normal = zeros(1, num_Normal);

% Find the PLV for all participants of both groups on both frequent and rare odors between
% the Fz and Cz channels using the function you implemented in section 3.4:
for p = 1:number_AD
    % Extract relevant data for the current AD patient
    person = AD(p);
    frequent = person.epoch(:,:,~logical(person.odor));
    rare = person.epoch(:,:,logical(person.odor));
    Fz_Cz_rare = squeeze(rare(2:3,:,:));
    Fz_Cz_frequent = squeeze(frequent(2:3,:,:));

    % Calculate PLV for rare trials
    Fz_mean = mean(Fz_Cz_rare(1,:,:),3);
    Cz_mean = mean(Fz_Cz_rare(2,:,:),3);
    plv_rare_AD(p) = PLV(Fz_mean,Cz_mean,sampling_rate,frequency_range);

    % Calculate PLV for frequent trials
    Fz_mean = mean(Fz_Cz_frequent(1,:,:),3);
    Cz_mean = mean(Fz_Cz_frequent(2,:,:),3);
    plv_frequent_AD(p) = PLV(Fz_mean,Cz_mean,sampling_rate,frequency_range);
end


% now for normal patients:
for p = 1:num_Normal
    % Extract relevant data for the current Normal patient
    person = normal(p);
    rare = person.epoch(:,:,logical(person.odor));
    frequent = person.epoch(:,:,~logical(person.odor));
    Fz_Cz_rare = squeeze(rare(2:3,:,:));
    Fz_Cz_frequent = squeeze(frequent(2:3,:,:));

    % Calculate PLV for rare trials
    Fz_mean = mean(Fz_Cz_rare(1,:,:),3);
    Cz_mean = mean(Fz_Cz_rare(2,:,:),3);
    plv_rare_Normal(p) = PLV(Fz_mean,Cz_mean,sampling_rate,frequency_range);

    % Calculate PLV for frequent trials
    Fz_mean = mean(Fz_Cz_frequent(1,:,:),3);
    Cz_mean = mean(Fz_Cz_frequent(2,:,:),3);
    plv_frequent_Normal(p) = PLV(Fz_mean,Cz_mean,sampling_rate,frequency_range);
end

%% 4.2

num_AD = numel(plv_frequent_AD);
num_Normal = numel(plv_frequent_Normal);

frequent_AD = [plv_frequent_AD, NaN(1, max_num - num_AD)].';
frequent_Normal = [plv_frequent_Normal, NaN(1, max_num - num_Normal)].';
rare_AD = [plv_rare_AD, NaN(1, max_num - num_AD)].';
rare_Normal = [plv_rare_Normal, NaN(1, max_num - num_Normal)].';

figure;
% boxplot
boxplot([frequent_AD, frequent_Normal, rare_AD, rare_Normal], ...
{'AD frequent', 'Normal frequent', 'AD rare', 'Normal rare'}, 'Symbol', 'o', 'Widths', 0.5);

title('PLV of AD and Normal Patients');
xlabel('Patient type');
ylabel('PLV');

figure;
pd = fitdist(frequent_AD, 'Normal');
x = linspace(0,1.3);
y = pdf(pd, x);
plot(x,y)

figure;
pd = fitdist(frequent_Normal, 'Normal');
x = linspace(0,1.3);
y = pdf(pd, x);
plot(x,y)

figure;
pd = fitdist(rare_AD, 'Normal');
x = linspace(0,1.3);
y = pdf(pd, x);
plot(x,y)

figure;
pd = fitdist(rare_Normal, 'Normal');
x = linspace(0,1.3);
y = pdf(pd, x);
plot(x,y)

%% 4.3
rare_normal_dist = fitdist(plv_rare_Normal', 'Normal');
rare_AD_dist = fitdist(plv_rare_AD', 'Normal');
frequent_normal_dist = fitdist(plv_frequent_Normal', 'Normal');
frequent_AD_dist = fitdist(plv_frequent_AD', 'Normal');

% t-test
[~, p_rare] = ttest2(plv_rare_Normal, plv_rare_AD, 'Vartype', 'unequal');
[~, p_freq] = ttest2(plv_frequent_Normal, plv_frequent_AD, 'Vartype', 'unequal');

distributions = {rare_normal_dist, rare_AD_dist, frequent_normal_dist, frequent_AD_dist};

for num = 1:numel(distributions)
    data = {plv_rare_Normal, plv_rare_AD, plv_frequent_Normal, plv_frequent_AD};
    data_index = data{num};
    a = linspace(min(data{num}), max(data{num}), 100);
    b = pdf(distributions{num}, a);
    figure;
    histogram(data_index, 6, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]);
    hold on;
    plot(a, b, 'LineWidth', 1);
    hold off;
    xlabel('PLV');
    ylabel('Prob');
    tit_name = {'rare-normal', 'rare-AD', 'Frequent-Normal', 'Frequent-AD'};
    title(tit_name{num});
end

%% 4.4
load('C:\Users\Soren\Desktop\Dataset\AD.mat');  
load('C:\Users\Soren\Desktop\Dataset\Normal.mat');  

patient_type = {normal, AD};

for i = 1:2
    type = patient_type{i};

    % Preallocate arrays to store phase differences
    n = numel(type) * size(type(1).epoch, 2);
    phaseDiff_freq = zeros(1, n);
    phaseDiff_rare = zeros(1, n);

    for i = 1:numel(type)
        
        person = type(i);
        rare = person.epoch(:, :, logical(person.odor));
        frequent = person.epoch(:, :, ~logical(person.odor));

        Fz_Cz_frequent = mean(squeeze(frequent(2:3, :, :)), 3);

        p = phase_diff(Fz_Cz_frequent(1,:), Fz_Cz_frequent(2,:));
        phaseDiff_freq((((i-1) * numel(p))+1):(((i-1) * numel(p))+numel(p))) = p;
        Fz_Cz_rare = mean(squeeze(rare(2:3, :, :)), 3);
        p = phase_diff(Fz_Cz_rare(1,:), Fz_Cz_rare(2,:));
        phaseDiff_rare((((i-1) * numel(p))+1):(((i-1) * numel(p))+numel(p))) = p;
    end

    figure;
    polarhistogram(phaseDiff_freq, 'Normalization', 'pdf');
    title('Phase-Diff');
    figure;
    polarhistogram(phaseDiff_rare, 'Normalization', 'pdf');
    title('Phase-Diff');
end























function plv = PLV(signals1, signals2, fs, frequencyrange)
    numSignals = size(signals1, 3);
    plv = zeros(1, numSignals);

    for i = 1:numSignals
        signal1 = signals1(:, :, i);
        signal2 = signals2(:, :, i);

        % Filtering the signal
        filteredSignal1 = bandpass(signal1, frequencyrange, fs);
        filteredSignal2 = bandpass(signal2, frequencyrange, fs);

        % Calculating phase difference
        analytic1 = hilbert(filteredSignal1);
        analytic2 = hilbert(filteredSignal2);
        phase1 = angle(analytic1);
        phase2 = angle(analytic2);
        phaseDiff = phase1 - phase2;

        % Calculate the PLV
        plv(i) = abs(mean(exp(1j * phaseDiff)));
    end
end

function phaseDiff = phase_diff(signal1, signal2)
    a1 = hilbert(signal1);
    a2 = hilbert(signal2);
    phase1 = angle(a1);
    phase2 = angle(a2);
    phaseDiff = phase1 - phase2;
end


