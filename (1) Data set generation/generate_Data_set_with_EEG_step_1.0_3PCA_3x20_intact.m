% Script for generating the X - Y data set,
% using all 1024 frequency components,
% in the data set the X not only contains the STFT spectrum, but also the EEG data of
% the first EEG window where the ACW is calculated from


number_of_components_to_use_each = 20; % !!!!!!!!!! subject to change
dataSetNum = 2; %%%!!!!! intact

step_size_between_EEG_windows = 1.0;

% measureNumList = [8,9,11,1];
measureNumList = [8];

outputFolder = 'Data_set_with_EEG_step_1.0_3PCA_3x20_intact';
mkdir(outputFolder)

% % load calculated STFT spectrum data to save some calculation time
% load(['CS545_data_set_20211208/x_matrix.mat']);
% x_matrix_20211208 = x_matrix;

% load PCA coeffs: 
% PCA_coeff, PCA_score, PCA_latent, PCA_tsquared, PCA_explained, PCA_mu
load('Spectrums - original and PCA/PCA_3xcoeffs_on_spectrogram_abs_r_i.mat')

dataFilesFolder = 'data/';
path_stimulus_resting = [dataFilesFolder, 'stimuli/', 'stim_122_123_pinkNoise.wav'];
path_stimulus_intact = [dataFilesFolder, 'stimuli/', 'stim_22_original.wav'];
path_stimulus_control = [dataFilesFolder, 'stimuli/', 'stim_23_phaseScrambledScaled.wav'];

load([dataFilesFolder, 'neural measure time series/', 'measure_ts_from_NMED_E_step_size_0_25.mat']); % !!!!!!
% shortTermAcousticTimeSeriesDir = [dataFilesFolder, '/NMts from stimuli/stim_23_phaseScrambledScaled.wav___shortTermFeatures.csv']; % !!!!!!

[data_stimulus_resting, ~] = audioread(path_stimulus_resting);
[data_stimulus_intact, Fs] = audioread(path_stimulus_intact);
[data_stimulus_control, ~] = audioread(path_stimulus_control);


% load the PCA coeffs of EEG data, 
% which contains eeg_PCA_coeff and eeg_PCA_explained
load('eeg_PCA_coeffs_whole.mat')

% Input EEG time series are stored in 'my_EEG_files'
path_data_file_122 = ['data/EEG/', ...
             'CleanEEG_stim122.mat'];
path_data_file_22 = ['data/EEG/', ...
             'CleanEEG_stim22.mat'];
path_data_file_123 = ['data/EEG/', ...
             'CleanEEG_stim123.mat'];      
path_data_file_23 = ['data/EEG/', ...
             'CleanEEG_stim23.mat'];


load(path_data_file_122);
load(path_data_file_22);
load(path_data_file_123);
load(path_data_file_23);

eeg_data_set{1} = eeg122;
eeg_data_set{2} = eeg22; % intact
eeg_data_set{3} = eeg123;
eeg_data_set{4} = eeg23; % control

clear eeg122 eeg22 eeg123 eeg23;



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
% get the list of neural measure names
subjectNumTotal = length(measureTimeSeriesFromEEG);
channelNumTotal = length(measureTimeSeriesFromEEG{1}{1});
table0 = measureTimeSeriesFromEEG{1}{1}{1}; % first subject; first channel; 
variableNameCellList = table0.Properties.VariableNames;
measureNameCellList = {};
for measureNum = 1:length(variableNameCellList)
    variableName = variableNameCellList{measureNum};
    measureName = variableName(7 : length(variableName)); % delete the "result" before each measure name
    measureNameCellList{measureNum} = measureName;
end
% measureNameCellList(2) = []; % for now, ignore the Spectral Entropy, since it cannot be calculated using the Cruncher


length_measureNumList = length(measureNumList);
% for number_measureNum = 1 : length_measureNumList
    
    number_measureNum = 1 ;
    
    measureNum = measureNumList(number_measureNum); % measure number
    fprintf('===========================================\n');
    fprintf('===========================================\n');
    fprintf('===========================================\n');
    fprintf('Start calculation for:\n');
    fprintf('%s, whose measureNum = %d\n', measureNameCellList{measureNum}, measureNum);
    % % measureNum = 8; % ACW
    % % measureNum = 9; % PLE
    % % measureNum = 11; % LZC
    % measureNum = 1; % MF
    % 
    % musicNum = 2;
    % subjectNum = 3;

    TN_musicNum = 1; % !!!!! total number of music
    TN_subjectNum = subjectNumTotal; % !!!!! total number of subjects
    
    resultX_matrices_by_loopNum = {};
    resultY_lists_by_loopNum = {};
    
    
    
    
%     subjectNum = 1;

for subjectNum = 1:TN_subjectNum
    
    fprintf('=======================================\n');
    fprintf('Processing subjectNum = %d\n\n', subjectNum);
    
    % deal with the EEG data
    eeg_data = eeg_data_set{dataSetNum};
    eeg_table_of_the_subject = eeg_data(:,:,subjectNum);
    
    
    
    stimulus_whole = data_stimulus_intact; %%%% !!!!!!!!!!!!!!!!!!!!!!!!!
    [stimulus_number_samples, ~] = size(stimulus_whole);


    h = height(measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{1});
    
    
    % count the number of samples from each subject
    num_samples_per_subject = 0;
    list_n = [];
    list_first_stimulus_window_num = [];
    list_last_music_window_num = [];
    n = 1; % data point number in EEG measure time series starts at #7 (i.e. from the 3rd second),
           % because we must leave some room for the former music windows
    first_stimulus_window_num = 1;
    stimulus_window_size = Fs * (3 + step_size_between_EEG_windows); %%%
    stimulus_window_hop_size = 11025; %%%
    last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
    
    
    % obtain the length of each column in x_matrix
    vec_stimulus_window = stimulus_whole(first_stimulus_window_num : last_music_window_num);
    spectrogram = [];
    % compute the spectrum of the stimulus window
    start_sample_num = 1;
    end_sample_num = start_sample_num + DFT_window_size - 1;
    while end_sample_num <= length(vec_stimulus_window)
        current_window = vec_stimulus_window(start_sample_num : end_sample_num);
        coef = B * current_window;
        spectrogram = [spectrogram, coef]; % using all 1024 frequency components

        start_sample_num = start_sample_num + DFT_hop_size;
        end_sample_num = start_sample_num + DFT_window_size - 1;
    end
    spectrogram_abs = abs(spectrogram);
    spectrogram_real = real(spectrogram);
    spectrogram_imag = imag(spectrogram);
    spectrogram_abs = spectrogram_abs'; % transpose for PCA
    spectrogram_real = spectrogram_real'; % transpose for PCA
    spectrogram_imag = spectrogram_imag'; % transpose for PCA
    
    spectrogram = spectrogram_abs;
    PCA_coeff = PCA_coeff_abs;
    % zero-mean normalization
    spectrogram_centered = spectrogram - mean(spectrogram, 1);
    % calculate the lower-dimensional representation of each sample
    spectrogram_PCA = spectrogram_centered * PCA_coeff(:, 1:number_of_components_to_use_each);
    % transpose back to the form we are familiar in the course
    % the number of rows equals to the number of PCA components used
    spectrogram_PCA_abs = spectrogram_PCA';
    
    spectrogram = spectrogram_real;
    PCA_coeff = PCA_coeff_real;
    % zero-mean normalization
    spectrogram_centered = spectrogram - mean(spectrogram, 1);
    % calculate the lower-dimensional representation of each sample
    spectrogram_PCA = spectrogram_centered * PCA_coeff(:, 1:number_of_components_to_use_each);
    % transpose back to the form we are familiar in the course
    % the number of rows equals to the number of PCA components used
    spectrogram_PCA_real = spectrogram_PCA';
    
    spectrogram = spectrogram_imag;
    PCA_coeff = PCA_coeff_imag;
    % zero-mean normalization
    spectrogram_centered = spectrogram - mean(spectrogram, 1);
    % calculate the lower-dimensional representation of each sample
    spectrogram_PCA = spectrogram_centered * PCA_coeff(:, 1:number_of_components_to_use_each);
    % transpose back to the form we are familiar in the course
    % the number of rows equals to the number of PCA components used
    spectrogram_PCA_imag = spectrogram_PCA';
    
    spectrogram_PCA = [spectrogram_PCA_abs; spectrogram_PCA_real; spectrogram_PCA_imag];
    
    eeg_starting_sample_num = floor(125 * (n - 1) / 4) + 1;
    eeg_ending_sample_num = eeg_starting_sample_num + 125 * 3 - 1; % 3-second EEG window
    matrix_eeg = eeg_table_of_the_subject(:, eeg_starting_sample_num:eeg_ending_sample_num);
    matrix_eeg = matrix_eeg';
    matrix_eeg_centered = matrix_eeg - mean(matrix_eeg, 1);
    % calculate the lower-dimensional representation of each sample
    matrix_eeg_PCA = matrix_eeg_centered * eeg_PCA_coeff(:, 1 : number_of_components_to_use_each * 3);
    matrix_eeg_PCA = matrix_eeg_PCA';
    
    matrix_all_features = [matrix_eeg_PCA, spectrogram_PCA];
    
    
    
    a_column_in_x_matrix = reshape(matrix_all_features, [], 1);
    
    num_rows_x_matrix = length(a_column_in_x_matrix);
    
    
    % now count the number of samples from each subject
    while ((n + 2) <= h) && (last_music_window_num <= stimulus_number_samples) %%%!!! stride is 0.5s
        num_samples_per_subject = num_samples_per_subject + 1;
        list_n = [list_n, n];
        list_first_stimulus_window_num = [list_first_stimulus_window_num, first_stimulus_window_num];
        list_last_music_window_num = [list_last_music_window_num, last_music_window_num];
        
        first_stimulus_window_num = first_stimulus_window_num + stimulus_window_hop_size * 2; %%%!!! stride is 0.5s
        last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
        
        n = n + 2; %%%!!! stride is 0.5s
    end
    
    fprintf('\nnum_samples_per_subject = %d\n\n', num_samples_per_subject);
    
    x_matrix = zeros(num_rows_x_matrix, num_samples_per_subject);
    y_list = zeros(num_samples_per_subject, 1);
    
    
    
    parfor sample_num = 1:num_samples_per_subject
        fprintf('Processing for sample #\t%d\n', sample_num);
        n = list_n(sample_num);
        first_stimulus_window_num = list_first_stimulus_window_num(sample_num);
        last_music_window_num = list_last_music_window_num(sample_num);
        
        % deal with X
        vec_stimulus_window = stimulus_whole(first_stimulus_window_num : last_music_window_num);

        spectrogram = [];
        % compute the spectrum of the stimulus window
        start_sample_num = 1;
        end_sample_num = start_sample_num + DFT_window_size - 1;
        while end_sample_num <= length(vec_stimulus_window)
            current_window = vec_stimulus_window(start_sample_num : end_sample_num);
            coef = B * current_window;
            spectrogram = [spectrogram, coef]; % using all 1024 frequency components

            start_sample_num = start_sample_num + DFT_hop_size;
            end_sample_num = start_sample_num + DFT_window_size - 1;
        end
        
        spectrogram_abs = abs(spectrogram);
        spectrogram_real = real(spectrogram);
        spectrogram_imag = imag(spectrogram);
        spectrogram_abs = spectrogram_abs'; % transpose for PCA
        spectrogram_real = spectrogram_real'; % transpose for PCA
        spectrogram_imag = spectrogram_imag'; % transpose for PCA

        spectrogram = spectrogram_abs;
        PCA_coeff = PCA_coeff_abs;
        % zero-mean normalization
        spectrogram_centered = spectrogram - mean(spectrogram, 1);
        % calculate the lower-dimensional representation of each sample
        spectrogram_PCA = spectrogram_centered * PCA_coeff(:, 1:number_of_components_to_use_each);
        % transpose back to the form we are familiar in the course
        % the number of rows equals to the number of PCA components used
        spectrogram_PCA_abs = spectrogram_PCA';

        spectrogram = spectrogram_real;
        PCA_coeff = PCA_coeff_real;
        % zero-mean normalization
        spectrogram_centered = spectrogram - mean(spectrogram, 1);
        % calculate the lower-dimensional representation of each sample
        spectrogram_PCA = spectrogram_centered * PCA_coeff(:, 1:number_of_components_to_use_each);
        % transpose back to the form we are familiar in the course
        % the number of rows equals to the number of PCA components used
        spectrogram_PCA_real = spectrogram_PCA';

        spectrogram = spectrogram_imag;
        PCA_coeff = PCA_coeff_imag;
        % zero-mean normalization
        spectrogram_centered = spectrogram - mean(spectrogram, 1);
        % calculate the lower-dimensional representation of each sample
        spectrogram_PCA = spectrogram_centered * PCA_coeff(:, 1:number_of_components_to_use_each);
        % transpose back to the form we are familiar in the course
        % the number of rows equals to the number of PCA components used
        spectrogram_PCA_imag = spectrogram_PCA';

        spectrogram_PCA = [spectrogram_PCA_abs; spectrogram_PCA_real; spectrogram_PCA_imag];
        
        % deal with the EEG
        eeg_starting_sample_num = floor(125 * (n - 1) / 4) + 1;
        eeg_ending_sample_num = eeg_starting_sample_num + 125 * 3 - 1; % 3-second EEG window
        matrix_eeg = eeg_table_of_the_subject(:, eeg_starting_sample_num:eeg_ending_sample_num);
        matrix_eeg = matrix_eeg';
        matrix_eeg_centered = matrix_eeg - mean(matrix_eeg, 1);
        % calculate the lower-dimensional representation of each sample
        matrix_eeg_PCA = matrix_eeg_centered * eeg_PCA_coeff(:, 1 : number_of_components_to_use_each * 3);
        matrix_eeg_PCA = matrix_eeg_PCA';

        matrix_all_features = [matrix_eeg_PCA, spectrogram_PCA];
        

        a_column_in_x_matrix = reshape(matrix_all_features, [], 1);
        
        
        x_matrix(:, sample_num) = a_column_in_x_matrix;


        % deal with Y
        SPC = 0;
        for electrodeNum = 1:125
            window_num_shift = step_size_between_EEG_windows * 4;
            w0 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n);
            w1 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n + window_num_shift); %%%!!! step size, not the stride
            percent_change = (w1 - w0) / w0;
            SPC = SPC + percent_change;
        end
        if SPC > 0
            y_list(sample_num) = 1;
        else
            y_list(sample_num) = 0;
        end
        
    end
    
    
    resultX_matrices_by_loopNum{subjectNum} = x_matrix;
    resultY_lists_by_loopNum{subjectNum} = y_list;
end

%     save([outputFolder, '/x_matrix.mat'], 'x_matrix', '-v7.3');
%     save('y_list', 'y_list', '-v7.3');


for subjectNum = 1:TN_subjectNum
    x_matrix = resultX_matrices_by_loopNum{subjectNum};
    y_list = resultY_lists_by_loopNum{subjectNum};
    path_X = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_X.mat'];
    path_Y = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_Y.mat'];
    save(path_X, 'x_matrix', '-v7.3');
    save(path_Y, 'y_list', '-v7.3');
end

% 
%     for subjectNum = 1:TN_subjectNum
% %         x_matrix = resultX_matrices_by_loopNum{subjectNum};
%         y_list = resultY_lists_by_loopNum{subjectNum};
% %         path_X = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_X.mat'];
%         path_Y = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_Y.mat'];
% %         save(path_X, 'x_matrix', '-v7.3');
%         save(path_Y, 'y_list', '-v7.3');
%     end

















