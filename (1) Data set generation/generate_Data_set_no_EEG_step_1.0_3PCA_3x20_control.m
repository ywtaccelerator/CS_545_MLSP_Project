% Script for generating the X - Y data set,
% using all 1024 frequency components,
% in the data set the X not only contains the STFT spectrum, but also the EEG data of
% the first EEG window where the ACW is calculated from


number_of_components_to_use_each = 20; % !!!!!!!!!! subject to change
dataSetNum = 4; %%%!!!!! control

step_size_between_EEG_windows = 1.0;

% measureNumList = [8,9,11,1];
measureNumList = [8];
% % measureNum = 8; % ACW
% % measureNum = 9; % PLE
% % measureNum = 11; % LZC
% measureNum = 1; % MF

outputFolder = 'Data_set_no_EEG_step_1.0_3PCA_3x20_control';
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

stimulus_whole = data_stimulus_control;  %%%% !!!!!!!!!!!!!!  control !!!
[stimulus_number_samples, ~] = size(stimulus_whole);




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
    
    
    subjectNum = 1;
    
    fprintf('Processing subjectNum = %d\n\n', subjectNum);
    
    

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
    a_column_in_x_matrix = reshape(spectrogram_PCA, [], 1);
    
    num_rows_x_matrix = length(a_column_in_x_matrix);
    
    
    % now count the number of samples from each subject
    while ((n + 1) <= h) && (last_music_window_num <= stimulus_number_samples) %%%!!! stride is 0.25s
        num_samples_per_subject = num_samples_per_subject + 1;
        list_n = [list_n, n];
        list_first_stimulus_window_num = [list_first_stimulus_window_num, first_stimulus_window_num];
        list_last_music_window_num = [list_last_music_window_num, last_music_window_num];
        
        first_stimulus_window_num = first_stimulus_window_num + stimulus_window_hop_size * 1; %%%!!! stride is 0.25s
        last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
        
        n = n + 1; %%%!!! stride is 0.25s
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
        a_column_in_x_matrix = reshape(spectrogram_PCA, [], 1);
        
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
    
    save([outputFolder, '/x_matrix.mat'], 'x_matrix', '-v7.3');
%     save('y_list', 'y_list', '-v7.3');


    
    
    
    parfor subjectNum = 2:TN_subjectNum %!! subject #1 is skipped since we have calculated above
            fprintf('Processing subjectNum = %d\n\n', subjectNum);

%             % read acoustic-feature matrix
%             shortTermFeatureTable = readtable(shortTermAcousticTimeSeriesDir);
%             
%             
%             table0VariableNameCellList = shortTermFeatureTable.Properties.VariableNames;
%             [numSamples, numFeatures] = size(shortTermFeatureTable);

%             stimulus_whole = data_stimulus_intact;
%             [stimulus_number_samples, ~] = size(stimulus_whole); %%%%

%             x_matrix = [];
            y_list = [];
            h = height(measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{1});
            
            
            n = 1; % data point number in EEG measure time series starts at #7 (i.e. from the 3rd second),
                   % because we must leave some room for the former music windows
            first_stimulus_window_num = 1;
%             stimulus_window_size = Fs * 4; %%%
%             stimulus_window_hop_size = 11025; %%%
            last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
            
            while ((n + 1) <= h) && (last_music_window_num <= stimulus_number_samples) %%%!!! stride is 0.25s
%                 % deal with X
%                 vec_stimulus_window = stimulus_whole(first_stimulus_window_num : last_music_window_num);
%                 
%                 spectrogram = [];
%                 % compute the spectrum of the stimulus window
%                 start_sample_num = 1;
%                 end_sample_num = start_sample_num + DFT_window_size - 1;
%                 while end_sample_num <= length(vec_stimulus_window)
%                     current_window = vec_stimulus_window(start_sample_num : end_sample_num);
%                     coef = B * current_window;
%                     spectrogram = [spectrogram, coef(1 : DFT_window_size / 2)];
% 
%                     start_sample_num = start_sample_num + DFT_hop_size;
%                     end_sample_num = start_sample_num + DFT_window_size - 1;
%                 end
%                 
%                 vec_spectrogram_real = reshape(real(spectrogram), [], 1);
%                 vec_spectrogram_imag = reshape(imag(spectrogram), [], 1);
%                 a_column_in_x_matrix = [vec_spectrogram_real; vec_spectrogram_imag];
%                 
%                 x_matrix = [x_matrix, a_column_in_x_matrix];
%                 
%                 
                first_stimulus_window_num = first_stimulus_window_num + stimulus_window_hop_size * 1; %%%!!! stride is 0.25s
                last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
                
                
                % deal with Y
                SPC = 0;
                for electrodeNum = 1:125
                    window_num_shift = step_size_between_EEG_windows * 4;
                    w0 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n);
                    w1 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n + window_num_shift); %%%!!!
                    percent_change = (w1 - w0) / w0;
                    SPC = SPC + percent_change;
                end

                if SPC > 0
                    y_list = [y_list; 1];
                else
                    y_list = [y_list; 0];
                end

                n = n + 1; %%%!!! stride is 0.25s
                
%                 if mod(n, 50) == 1
%                     fprintf('%d out of %d was processed\n', floor(n/2) + 1, floor(h/2) + 1);
%                 end
            end
            

%             resultX_matrices_by_loopNum{subjectNum} = x_matrix;
            resultY_lists_by_loopNum{subjectNum} = y_list;

%                 resultX_matrix = [resultX_matrix; x_matrix];
%                 resultY_list = [resultY_list; y_list];
%             end
%         end
    end
    
%     for subjectNum = 1:TN_subjectNum
%         x_matrix = resultX_matrices_by_loopNum{subjectNum};
%         y_list = resultY_lists_by_loopNum{subjectNum};
%         path_X = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_X.mat'];
%         path_Y = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_Y.mat'];
%         save(path_X, 'x_matrix', '-v7.3');
%         save(path_Y, 'y_list', '-v7.3');
%     end
    
    for subjectNum = 1:TN_subjectNum
%         x_matrix = resultX_matrices_by_loopNum{subjectNum};
        y_list = resultY_lists_by_loopNum{subjectNum};
%         path_X = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_X.mat'];
        path_Y = [outputFolder, '/sj_', num2str(subjectNum, '%02d'), '_Y.mat'];
%         save(path_X, 'x_matrix', '-v7.3');
        save(path_Y, 'y_list', '-v7.3');
    end

















