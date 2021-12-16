% do PCA on the STFT spectrum of all three stimuli:
% calculate the STFT spectrum of all three stimuli
% do PCA on the spectrogram of a whole experiment, and save all PCA coeffs

outputFolder = 'Spectrums - original and PCA';
mkdir(outputFolder)



dataFilesFolder = 'data/';
path_stimulus_resting = [dataFilesFolder, 'stimuli/', 'stim_122_123_pinkNoise.wav'];
path_stimulus_intact = [dataFilesFolder, 'stimuli/', 'stim_22_original.wav'];
path_stimulus_control = [dataFilesFolder, 'stimuli/', 'stim_23_phaseScrambledScaled.wav'];

load([dataFilesFolder, 'neural measure time series/', 'measure_ts_from_NMED_E_step_size_0_25.mat']); % !!!!!!
% shortTermAcousticTimeSeriesDir = [dataFilesFolder, '/NMts from stimuli/stim_23_phaseScrambledScaled.wav___shortTermFeatures.csv']; % !!!!!!

[data_stimulus_resting, ~] = audioread(path_stimulus_resting);
[data_stimulus_intact, Fs] = audioread(path_stimulus_intact);
[data_stimulus_control, ~] = audioread(path_stimulus_control);



% prepare for DFT
DFT_window_size = 1024;
N = DFT_window_size;
DFT_hop_size = 512;
F = [];
for j = 1:N
    for k  = 1:N
        F(j,k) = (1/sqrt(N)) * exp(1i*j*k*2*pi/N);
    end
end
w = 0.5.*(1 - cos(2*pi.*((1:N)-1)./(N-1))); % hamming window
H = diag(w);

B = F * H;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%!!!!! process for intact
stimulus_whole = data_stimulus_intact; 


spectrogram = [];
% compute the spectrum of the stimulus window
start_sample_num = 1;
end_sample_num = start_sample_num + DFT_window_size - 1;
stimulus_length = length(stimulus_whole);
percent_processed_print = 0;
while end_sample_num <= stimulus_length
    percent_processed = floor(start_sample_num / stimulus_length * 100);
    if (percent_processed_print == percent_processed) && mod(percent_processed, 5) == 0
        fprintf('Processed %d percent\n', percent_processed);
        percent_processed_print = percent_processed_print + 5;
    end
    current_window = stimulus_whole(start_sample_num : end_sample_num);
    coef = B * current_window;
    spectrogram = [spectrogram, coef]; % using all 1024 frequency components

    start_sample_num = start_sample_num + DFT_hop_size;
    end_sample_num = start_sample_num + DFT_window_size - 1;
end

spectrogram_intact = spectrogram;
save('spectrogram_complex_intact.mat', 'spectrogram_intact', '-v7.3');



%%%!!!!! process for control
stimulus_whole = data_stimulus_control; 


spectrogram = [];
% compute the spectrum of the stimulus window
start_sample_num = 1;
end_sample_num = start_sample_num + DFT_window_size - 1;
stimulus_length = length(stimulus_whole);
percent_processed_print = 0;
while end_sample_num <= stimulus_length
    percent_processed = floor(start_sample_num / stimulus_length * 100);
    if (percent_processed_print == percent_processed) && mod(percent_processed, 5) == 0
        fprintf('Processed %d percent\n', percent_processed);
        percent_processed_print = percent_processed_print + 5;
    end
    current_window = stimulus_whole(start_sample_num : end_sample_num);
    coef = B * current_window;
    spectrogram = [spectrogram, coef]; % using all 1024 frequency components

    start_sample_num = start_sample_num + DFT_hop_size;
    end_sample_num = start_sample_num + DFT_window_size - 1;
end

spectrogram_control = spectrogram;
save('spectrogram_complex_control.mat', 'spectrogram_control', '-v7.3');



%%%!!!!! process for resting
stimulus_whole = data_stimulus_resting; 


spectrogram = [];
% compute the spectrum of the stimulus window
start_sample_num = 1;
end_sample_num = start_sample_num + DFT_window_size - 1;
stimulus_length = length(stimulus_whole);
percent_processed_print = 0;
while end_sample_num <= stimulus_length
    percent_processed = floor(start_sample_num / stimulus_length * 100);
    if (percent_processed_print == percent_processed) && mod(percent_processed, 5) == 0
        fprintf('Processed %d percent\n', percent_processed);
        percent_processed_print = percent_processed_print + 5;
    end
    current_window = stimulus_whole(start_sample_num : end_sample_num);
    coef = B * current_window;
    spectrogram = [spectrogram, coef]; % using all 1024 frequency components

    start_sample_num = start_sample_num + DFT_hop_size;
    end_sample_num = start_sample_num + DFT_window_size - 1;
end

spectrogram_resting = spectrogram;
save('spectrogram_complex_resting.mat', 'spectrogram_resting', '-v7.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% do PCA on the spectrogram of a whole experiment

load('Spectrums - original and PCA/spectrogram_complex_intact.mat')
load('Spectrums - original and PCA/spectrogram_complex_control.mat')
load('Spectrums - original and PCA/spectrogram_complex_resting.mat')

% for intact
spectrogram = spectrogram_intact;
spectrogram_whole_real = real(spectrogram);
spectrogram_whole_imag = imag(spectrogram);
spectrogram_whole_r_i = [spectrogram_whole_real; spectrogram_whole_imag];
spectrogram_whole_r_i_T_intact = spectrogram_whole_r_i'; % transpose for PCA

% for control
spectrogram = spectrogram_control;
spectrogram_whole_real = real(spectrogram);
spectrogram_whole_imag = imag(spectrogram);
spectrogram_whole_r_i = [spectrogram_whole_real; spectrogram_whole_imag];
spectrogram_whole_r_i_T_control = spectrogram_whole_r_i'; % transpose for PCA

% for resting
spectrogram = spectrogram_resting;
spectrogram_whole_real = real(spectrogram);
spectrogram_whole_imag = imag(spectrogram);
spectrogram_whole_r_i = [spectrogram_whole_real; spectrogram_whole_imag];
spectrogram_whole_r_i_T_resting = spectrogram_whole_r_i'; % transpose for PCA


% spectrogram of a whole experiment
spectrogram_whole_r_i_T = [
    spectrogram_whole_r_i_T_resting;
    spectrogram_whole_r_i_T_intact;
    spectrogram_whole_r_i_T_resting;
    spectrogram_whole_r_i_T_control];



[PCA_coeff, PCA_score, PCA_latent, PCA_tsquared, PCA_explained, PCA_mu] = pca(spectrogram_whole_r_i_T);

save('PCA_coeffs_on_spectrogram_whole.mat', 'PCA_coeff', 'PCA_score', 'PCA_latent', 'PCA_tsquared', 'PCA_explained', 'PCA_mu', '-v7.3');



% ------------------------------------------------------------------------------------------

% number_of_components_to_use = 


percent_explained_threshold = 95;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 59 when percent_explained_threshold = 95
% let us pick number_of_components_to_use = 60
% ------------------------------------------------------------------------------------------
percent_explained_threshold = 98;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 94 when percent_explained_threshold = 98
% let us pick number_of_components_to_use = 100
% ------------------------------------------------------------------------------------------
percent_explained_threshold = 99;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 143 when percent_explained_threshold = 99
% let us pick number_of_components_to_use = 150
% ------------------------------------------------------------------------------------------
percent_explained_threshold = 90;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 40 when percent_explained_threshold = 90
% let us pick number_of_components_to_use = 40
% ------------------------------------------------------------------------------------------
percent_explained_threshold = 85;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 31 when percent_explained_threshold = 85
% let us pick number_of_components_to_use = 30






pareto(latent1);%调用matla画图 pareto仅绘制累积分布的前95%，因此y中的部分元素并未显示
xlabel('Principal Component');
ylabel('Variance Explained (%)');
% 图中的线表示的累积变量解释程度
print(gcf,'-dpng','PCA.png');
sum_explained = sum(PCA_explained(1 : number_of_components_to_use));





[data_PCA, coeff, sum_explained, number_of_components_to_use]=pca_demo_2(data)
%用percent_threshold决定保留xx%的贡献率
percent_explained_threshold = 95;   %百分比阈值，用于决定保留的主成分个数；
% data=zscore(data);  %归一化数据
[coeff, score, latent, tsquared, explained, mu] = pca(data);
latent1=100*latent/sum(latent);%将latent特征值总和统一为100，便于观察贡献率
A=length(latent1);
percents=0;                          %累积百分比
for number_of_components_to_use = 1:A
    percents=percents+latent1(number_of_components_to_use);
    if percents > percent_explained_threshold
        break;
    end
end
% data= bsxfun(@minus,data,mean(data,1));
data_PCA = data * coeff(:, 1:number_of_components_to_use);
pareto(latent1);%调用matla画图 pareto仅绘制累积分布的前95%，因此y中的部分元素并未显示
xlabel('Principal Component');
ylabel('Variance Explained (%)');
% 图中的线表示的累积变量解释程度
print(gcf,'-dpng','PCA.png');
sum_explained=sum(explained(1:number_of_components_to_use));









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% do PCA on the magnitude/real/imaginary spectrogram of a whole experiment

load('Spectrums - original and PCA/spectrogram_complex_intact.mat')
load('Spectrums - original and PCA/spectrogram_complex_control.mat')
load('Spectrums - original and PCA/spectrogram_complex_resting.mat')

% for intact
spectrogram = spectrogram_intact;
spectrogram_whole_real = real(spectrogram);
spectrogram_whole_imag = imag(spectrogram);
spectrogram_whole_r_i = [spectrogram_whole_real; spectrogram_whole_imag];
spectrogram_whole_r_i_T_intact = spectrogram_whole_r_i'; % transpose for PCA

% for control
spectrogram = spectrogram_control;
spectrogram_whole_real = real(spectrogram);
spectrogram_whole_imag = imag(spectrogram);
spectrogram_whole_r_i = [spectrogram_whole_real; spectrogram_whole_imag];
spectrogram_whole_r_i_T_control = spectrogram_whole_r_i'; % transpose for PCA

% for resting
spectrogram = spectrogram_resting;
spectrogram_whole_real = real(spectrogram);
spectrogram_whole_imag = imag(spectrogram);
spectrogram_whole_r_i = [spectrogram_whole_real; spectrogram_whole_imag];
spectrogram_whole_r_i_T_resting = spectrogram_whole_r_i'; % transpose for PCA


% spectrogram of a whole experiment
spectrogram_whole_r_i_T = [
    spectrogram_whole_r_i_T_resting;
    spectrogram_whole_r_i_T_intact;
    spectrogram_whole_r_i_T_resting;
    spectrogram_whole_r_i_T_control];







spectrogram_whole_complex = [
    spectrogram_resting';
    spectrogram_intact';
    spectrogram_resting';
    spectrogram_control'];

spectrogram_whole_abs = abs(spectrogram_whole_complex);
spectrogram_whole_real = real(spectrogram_whole_complex);
spectrogram_whole_imag = imag(spectrogram_whole_complex);



[PCA_coeff, PCA_score, PCA_latent, PCA_tsquared, PCA_explained, PCA_mu] = pca(spectrogram_whole_abs);

PCA_coeff_abs = PCA_coeff;
PCA_explained_abs = PCA_explained;

[PCA_coeff, PCA_score, PCA_latent, PCA_tsquared, PCA_explained, PCA_mu] = pca(spectrogram_whole_real);

PCA_coeff_real = PCA_coeff;
PCA_explained_real = PCA_explained;

[PCA_coeff, PCA_score, PCA_latent, PCA_tsquared, PCA_explained, PCA_mu] = pca(spectrogram_whole_imag);

PCA_coeff_imag = PCA_coeff;
PCA_explained_imag = PCA_explained;


matrix_PCA_explained = [PCA_explained_abs, PCA_explained_real, PCA_explained_imag];
save('matrix_PCA_explained.mat', 'matrix_PCA_explained', '-v7.3');

save('PCA_3xcoeffs_on_spectrogram_abs_r_i.mat', ...
    'PCA_coeff_abs', 'PCA_explained_abs', ...
    'PCA_coeff_real', 'PCA_explained_real', ...
    'PCA_coeff_imag', 'PCA_explained_imag', '-v7.3');

% ------------------------------------------------------------------------------------------

% number_of_components_to_use = 

PCA_explained = PCA_explained_abs;

percent_explained_threshold = 95;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 59 
% when percent_explained_threshold = 95
% let us pick number_of_components_to_use = 60

sum_explained=sum(PCA_explained(1:number_of_components_to_use));
% to explain 95%, we need abs-real-imag = 25-30-29
% so we pick 30 from each of them

% ------------------------------------------------------------------------------------------
PCA_explained = PCA_explained_imag;

percent_explained_threshold = 90;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end

% to explain 90%, we need abs-real-imag = 16-20-20
% so we pick 20 from each of them

% ------------------------------------------------------------------------------------------
PCA_explained = PCA_explained_imag;

percent_explained_threshold = 85;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 94 
% when percent_explained_threshold = 98
% let us pick number_of_components_to_use = 100
sum_explained=sum(PCA_explained_real(1:15));
% to explain 85%, we need abs-real-imag = 11-16-15
% so we pick 15 from each of them

% ------------------------------------------------------------------------------------------
PCA_explained = PCA_explained_imag;

percent_explained_threshold = 80;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 94 
% when percent_explained_threshold = 98
% let us pick number_of_components_to_use = 100
sum_explained=sum(PCA_explained_real(1:15));
% to explain 85%, we need abs-real-imag = 8-13-12
% so we pick 12 from each of them

% ------------------------------------------------------------------------------------------
PCA_explained = PCA_explained_real;

percent_explained_threshold = 98;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
% we can get number_of_components_to_use = 59 
% when percent_explained_threshold = 95
% let us pick number_of_components_to_use = 60

sum_explained=sum(PCA_explained(1:number_of_components_to_use));
% to explain 98%, we need abs-real-imag = 41-47-47
% so we pick 50 from each of them





[PCA_coeff, PCA_score, PCA_latent, PCA_tsquared, PCA_explained, PCA_mu] = pca(spectrogram_whole_r_i_T);

save('PCA_coeffs_on_spectrogram_whole.mat', 'PCA_coeff', 'PCA_score', 'PCA_latent', 'PCA_tsquared', 'PCA_explained', 'PCA_mu', '-v7.3');











% for EEG:

eeg_table_all = [];
for subjectNum = 1:TN_subjectNum
    for data_set_num = 1:4
        eeg_data = eeg_data_set{data_set_num};
        eeg_table_of_the_subject = eeg_data(:,:,subjectNum);
        eeg_table_of_the_subject = eeg_table_of_the_subject';
        eeg_table_all = [eeg_table_all; eeg_table_of_the_subject];
    end
end

[eeg_PCA_coeff, ~, ~, ~, eeg_PCA_explained, ~] = pca(eeg_table_all);

save('eeg_PCA_coeffs_whole.mat', 'eeg_PCA_coeff', 'eeg_PCA_explained', '-v7.3');


% 90-95-98:
% 42-65-91
% ------------------------------------------------------------------------------------------
PCA_explained = eeg_PCA_explained;

percent_explained_threshold = 95;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
number_of_components_to_use
% we can get number_of_components_to_use = 65 
% when percent_explained_threshold = 95
% let us pick number_of_components_to_use = 100
sum_explained=sum(PCA_explained_real(1:15));


% ------------------------------------------------------------------------------------------
PCA_explained = eeg_PCA_explained;

percent_explained_threshold = 98;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
number_of_components_to_use
% we can get number_of_components_to_use = 91
% when percent_explained_threshold = 98

% ------------------------------------------------------------------------------------------
PCA_explained = eeg_PCA_explained;

percent_explained_threshold = 90;
% obtain the number_of_components_to_use
num_features = length(PCA_explained);
accumulated_percent = 0;
for number_of_components_to_use = 1:num_features
    accumulated_percent = accumulated_percent + PCA_explained(number_of_components_to_use);
    if accumulated_percent > percent_explained_threshold
        break;
    end
end
number_of_components_to_use
% we can get number_of_components_to_use = 42
% when percent_explained_threshold = 90






































































































































































































