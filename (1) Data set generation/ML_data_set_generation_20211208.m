% Script for generating the X - Y data set,
% using all 1024 frequency components,
% in the data set the X not only contains the STFT spectrum, but also the EEG data of
% the first EEG window where the ACW is calculated from


dataSetNum = 2; %%%!!!!! intact


% measureNumList = [8,9,11,1];
measureNumList = [8];

outputFolder = 'CS545_data_set_20211208';
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
for number_measureNum = 1 : length_measureNumList
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
    
    stimulus_whole = data_stimulus_intact; %%%%
    [stimulus_number_samples, ~] = size(stimulus_whole);

    x_matrix = [];
    y_list = [];
    h = height(measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{1});


    n = 1; % data point number in EEG measure time series starts at #7 (i.e. from the 3rd second),
           % because we must leave some room for the former music windows
    first_stimulus_window_num = 1;
    stimulus_window_size = Fs * 3.5; %%%
    stimulus_window_hop_size = 11025; %%%
    last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;

    while ((n + 2) <= h) && (last_music_window_num <= stimulus_number_samples) %%%!!!
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

        vec_spectrogram_real = reshape(real(spectrogram), [], 1);
        vec_spectrogram_imag = reshape(imag(spectrogram), [], 1);
        a_column_in_x_matrix = [vec_spectrogram_real; vec_spectrogram_imag];

        x_matrix = [x_matrix, a_column_in_x_matrix];


        first_stimulus_window_num = first_stimulus_window_num + stimulus_window_hop_size * 2; %%%!!! stride is 0.5s
        last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;


        % deal with Y
        SPC = 0;
        for electrodeNum = 1:125
            w0 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n);
            w1 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n + 2); %%%!!!
            percent_change = (w1 - w0) / w0;
            SPC = SPC + percent_change;
        end

        if SPC > 0
            y_list = [y_list; 1];
        else
            y_list = [y_list; 0];
        end

        n = n + 2; %%%!!! stride is 0.5s

        if mod(n, 50) == 1
            fprintf('%d out of %d was processed\n', floor(n/2) + 1, floor(h/2) + 1);
        end
    end
    
    resultX_matrices_by_loopNum{subjectNum} = x_matrix;
    resultY_lists_by_loopNum{subjectNum} = y_list;
    
    save('x_matrix', 'x_matrix', '-v7.3');
    save('y_list', 'y_list', '-v7.3');
    
    
    
    for subjectNum = 2:TN_subjectNum %!! subject #1 is skipped since we have calculated above
            fprintf('Processing subjectNum = %d\n\n', subjectNum);

%             % read acoustic-feature matrix
%             shortTermFeatureTable = readtable(shortTermAcousticTimeSeriesDir);
%             
%             
%             table0VariableNameCellList = shortTermFeatureTable.Properties.VariableNames;
%             [numSamples, numFeatures] = size(shortTermFeatureTable);

            stimulus_whole = data_stimulus_intact;
            [stimulus_number_samples, ~] = size(stimulus_whole); %%%%

%             x_matrix = [];
            y_list = [];
            h = height(measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{1});
            
            
            n = 1; % data point number in EEG measure time series starts at #7 (i.e. from the 3rd second),
                   % because we must leave some room for the former music windows
            first_stimulus_window_num = 1;
            stimulus_window_size = Fs * 4; %%%
            stimulus_window_hop_size = 11025; %%%
            last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
            
            while ((n + 2) <= h) && (last_music_window_num <= stimulus_number_samples) %%%!!!
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
%                 first_stimulus_window_num = first_stimulus_window_num + stimulus_window_hop_size * 2; %%%!!! stride is 0.5s
%                 last_music_window_num = first_stimulus_window_num + stimulus_window_size - 1;
                
                
                % deal with Y
                SPC = 0;
                for electrodeNum = 1:125
                    w0 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n);
                    w1 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n + 2); %%%!!!
                    percent_change = (w1 - w0) / w0;
                    SPC = SPC + percent_change;
                end

                if SPC > 0
                    y_list = [y_list; 1];
                else
                    y_list = [y_list; 0];
                end

                n = n + 2; %%%!!! stride is 0.5s
                
                if mod(n, 50) == 1
                    fprintf('%d out of %d was processed\n', floor(n/2) + 1, floor(h/2) + 1);
                end
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
    
    
    % (optional) shortern the x_matrix
    [~, width_x_matrix] = size(x_matrix);
    width_x_matrix_reduced = [];
    for column_num = 1:width_x_matrix
        x = x_matrix(:, column_num);
        x_r = reshape(x, 512, []);
        x_r_reduced = x_r(1:160, :);
        vec_x_r_reduced = reshape(x_r_reduced, [], 1);
        width_x_matrix_reduced = [width_x_matrix_reduced, vec_x_r_reduced];
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    data_set.training_set_X = training_set_X_cells;
    data_set.training_set_Y = training_set_Y_cells;
    data_set.test_set_X = test_set_X_cells;
    data_set.test_set_Y = test_set_Y_cells;

    outputPath = [outputFolder, '/EEG_SAPC_data_set.mat'];
    save(outputPath,'data_set', '-v7.3');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    subFolderPath1 = [path1, path2, path3, path4, '_', path_measureName, '/'];
    mkdir(subFolderPath1);
    fileName1 = ['vA_Ycr2_pct', path2, path3, path4, '_', path_measureName];

    infotxtPath = [subFolderPath1, fileName1, '_info.txt'];
    infotxt_fid = fopen(infotxtPath,'w');
    fprintf(infotxt_fid, 'Experiment Parameter:\n%s\n\n', fileName1);

    count0_total = 0;
    count1_total = 0;

    for subjectNum = 1 : TN_subjectNum
        fprintf(infotxt_fid, 'Subject # %d:\n', subjectNum);
        resultX_matrix = [];
        resultY_list = [];
        for musicNum = 1 : TN_musicNum
%             fprintf("%d\n", look_up_table{subjectNum}{musicNum});
            loopNum = look_up_table{subjectNum}{musicNum};
            X_matrix = resultX_matrices_by_loopNum{loopNum};
            Y_list = resultY_lists_by_loopNum{loopNum};
            resultX_matrix = [resultX_matrix; X_matrix];
            resultY_list = [resultY_list; Y_list];
        end

        count0 = 0;
        count1 = 0;
        for i = 1 : length(resultY_list)
            if resultY_list(i) == 1
                count1 = count1 + 1;
            else
                count0 = count0 + 1;
            end
        end

        count0_total = count0_total + count0;
        count1_total = count1_total + count1;

        % output the actual number of 1??s and the actual number of 0??s
        fprintf(infotxt_fid, 'The actual number of Y=1''s to the actual number of Y=0''s:\n');
        fprintf(infotxt_fid, '[%d:%d]\n', count1, count0);
        % output the ratio of the number of Y=1's to the number of Y=0's
        fprintf(infotxt_fid, 'The ratio of the number of Y=1''s to the number of Y=0''s:\n');
        fprintf(infotxt_fid, '%.8f\n\n', (count1 / count0));


        % save("EEGmeasures_longTermFeature_for_mlp/resultX_table.mat", "resultX_table");

        % writetable(resultX_table, ...
        %     "EEGmeasures_longTermFeature_for_mlp/vA_Ycr2_pct25_ep_2_3_MF_resultX_table.csv");
        % save("EEGmeasures_longTermFeature_for_mlp/vA_Ycr2_pct25_ep_2_3_MF_resultY_list.mat", "resultY_list");

        subFolderPath2 = [subFolderPath1, 'sj', num2str(subjectNum, '%02d'), '/myData20200815/'];
        mkdir(subFolderPath2);
        path_X = [subFolderPath2, fileName1, '_resultX_matrix.mat'];
        path_Y = [subFolderPath2, fileName1, '_resultY_list.mat'];

%         writetable(resultX_table, path_X);
        save(path_X, 'resultX_matrix', '-v7.3');
        save(path_Y, 'resultY_list');

    end

    % output the statistical information:
    fprintf(infotxt_fid, '============================================================\n');
    fprintf(infotxt_fid, 'In total, under this experiment parameter:\n');
    fprintf(infotxt_fid, 'The actual number of Y=1''s to the actual number of Y=0''s:\n');
    fprintf(infotxt_fid, '[%d:%d]\n', count1_total, count0_total);
    fprintf(infotxt_fid, 'The ratio of the number of Y=1''s to the number of Y=0's:\n');
    fprintf(infotxt_fid, '%.8f\n\n', (count1_total / count0_total));

    fclose(infotxt_fid);
end


































































































































































