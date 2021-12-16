
dataSetNum = 2; %%%!!!!! intact


% measureNumList = [8,9,11,1];
measureNumList = [8];

outputFolder = 'CS545_data_set_20211120';
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

    parfor subjectNum = 1:TN_subjectNum
            fprintf('Processing subjectNum = %d\n\n', subjectNum);

%             % read acoustic-feature matrix
%             shortTermFeatureTable = readtable(shortTermAcousticTimeSeriesDir);
%             
%             
%             table0VariableNameCellList = shortTermFeatureTable.Properties.VariableNames;
%             [numSamples, numFeatures] = size(shortTermFeatureTable);

            [stimulus_number_samples, ~] = size(data_stimulus_intact); %%%%

            x_matrix = [];
            y_list = [];
            h = height(measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{1});
%             h_music = height(shortTermFeatureTable);
            
            
            
            n = 1; % data point number in EEG measure time series starts at #7 (i.e. from the 3rd second),
                   % because we must leave some room for the former music windows
            while ((n + 4) <= h)
                
                
                
                
                

                % deal with X
                first_music_window_number = 20 * (n - 1) + 1;
                last_music_window_number = 20 * (n - 1) + 239;
                a_row_of_resultX_matrix = [];
                music_part = [];
                doesHaveNaN = false;
                for music_window_number = first_music_window_number:last_music_window_number
                    a_row_of_shortTermFeatureTable = shortTermFeatureTable{music_window_number,:};
                    numNaN = numel(find(isnan(a_row_of_shortTermFeatureTable)));
                    if numNaN > 0
                        % if a_row_of_resultX_matrix has a NaN, we
                        % simply do not use this row of resultX_matrix
                        doesHaveNaN = true;
                        break;
                    end

                    music_part = [music_part, a_row_of_shortTermFeatureTable'];
%                         % show the values of X for the current
%                         % W(old)-W(new) pair
%                         disp(reshape(a_row_of_resultX_matrix, 5, 13)')
%                         disp(reshape(a_row_of_resultX_matrix, 5, 13))
                end

                if doesHaveNaN == true
                    n = n + 1;
                    continue;
                else
                    % the time stamp where the window of a neural measure
                    % starts
                    time_stamp_of_n = (n - 1) / 4;
                    
                    % the number of the sample in raw EEG that starts from
                    % the time stamp above
                    starting_sample_num_in_raw_EEG = round(time_stamp_of_n * 125 + 1);
                    ending_sample_num_in_raw_EEG = starting_sample_num_in_raw_EEG + 375 - 1;
                    
                    eeg_raw_data = eeg_raw_data_original(:,:,subjectNum);
                    
                    [sr_raw_EEG , h_raw_EEG] = size(eeg_raw_data);
                    
                    if ending_sample_num_in_raw_EEG > h_raw_EEG
                        n = n + 1;
                        continue;
                    else
                        % enlarge the music part of the input X
                        doubled_music_part = [];
                        [row_num_total, column_num_total] = size(music_part);
                        for row_num = 1:row_num_total
                            a_row_in_doubled_music_part = [];
                            for column_num = 1:column_num_total
                                element = music_part(row_num, column_num);
                                group = [element, element];
                                a_row_in_doubled_music_part = [a_row_in_doubled_music_part, group];

                            end
                            doubled_music_part = [doubled_music_part; a_row_in_doubled_music_part];
                            doubled_music_part = [doubled_music_part; a_row_in_doubled_music_part];
                        end
                        
                        [height_doubled_music_part, width_doubled_music_part] = size(doubled_music_part);
                        
                        raw_EEG_part = eeg_raw_data(:, starting_sample_num_in_raw_EEG:ending_sample_num_in_raw_EEG);
                        [height_raw_EEG_part, width_raw_EEG_part] = size(raw_EEG_part);
                        
                        zero_padding = zeros(125, width_doubled_music_part - width_raw_EEG_part);
                        padded_raw_EEG_part = [raw_EEG_part, zero_padding];
                        
                        final_input_X = [padded_raw_EEG_part; doubled_music_part];
                        
                        a_row_of_resultX_matrix = reshape(final_input_X', 1, []);
                    end
                    
                    
                    
                    x_matrix = [x_matrix; a_row_of_resultX_matrix];
                    % x_matrix stores all of the X samples of a subject
                    % under a music, each sample is a row of x_matrix
                end


                % deal with Y
            %     w0w1s = [];
                SAPC = 0;
                for electrodeNum = 1:125
                    w0 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n);
                    w1 = measureTimeSeriesFromEEG{subjectNum}{dataSetNum}{electrodeNum}{:, measureNum}(n + 4);
            %         w0w1s(electrodeNum,:) = [w0, w1];
                    percent_change = (w1 - w0) / w0;
                    SAPC = SAPC + percent_change;
                end

                if SAPC > 0
                    y_list = [y_list; 1];
                else
                    y_list = [y_list; 0];
                end

                n = n + 1;
            end


            resultX_matrices_by_loopNum{subjectNum} = x_matrix;
            resultY_lists_by_loopNum{subjectNum} = y_list;

%                 resultX_matrix = [resultX_matrix; x_matrix];
%                 resultY_list = [resultY_list; y_list];
%             end
%         end
    end

    
    % generate the path
    if measureNum == 8
        path_measureName = 'ACW';
    else
        if measureNum == 9
            path_measureName = 'PLE';
        else
            if measureNum == 11
                path_measureName = 'LZC';
            else
                if measureNum == 1
                    path_measureName = 'MF';
                end
            end
        end
    end
    
    
    
    
    sp_tr_1_start = 1;
    sp_tr_1_end = 109;
    sp_tr_2_start = 181;
    sp_tr_2_end = 709;
    sp_tr_3_start = 781;
    sp_tr_3_end = 1429;
    sp_tr_4_start = 1501;
    
    sp_te_1_start = 121;
    sp_te_1_end = 169;
    sp_te_2_start = 721;
    sp_te_2_end = 769;
    sp_te_3_start = 1441;
    sp_te_3_end = 1489;
    
    training_set_X_cells = [];
    training_set_Y_cells = [];
    test_set_X_cells = [];
    test_set_Y_cells = [];
    
    for subjectNum = 1 : TN_subjectNum
        fprintf('Processing subject #%d:\n', subjectNum);
        
        for musicNum = 1 : TN_musicNum
            fprintf('Processing music #%d:\n', musicNum);
%             fprintf("%d\n", look_up_table{subjectNum}{musicNum});
            loopNum = subjectNum;
            X_matrix = resultX_matrices_by_loopNum{loopNum};
            Y_list = resultY_lists_by_loopNum{loopNum};
            
            [h_X_matrix, w_X_matrix] = size(X_matrix);
            
            test_X_p1 = X_matrix(sp_te_1_start: sp_te_1_end, :);
            test_X_p2 = X_matrix(sp_te_2_start: sp_te_2_end, :);
            test_X_p3 = X_matrix(sp_te_3_start: sp_te_3_end, :);
            test_X = [test_X_p1; test_X_p2; test_X_p3];
            
            test_Y_p1 = Y_list(sp_te_1_start: sp_te_1_end, :);
            test_Y_p2 = Y_list(sp_te_2_start: sp_te_2_end, :);
            test_Y_p3 = Y_list(sp_te_3_start: sp_te_3_end, :);
            test_Y = [test_Y_p1; test_Y_p2; test_Y_p3];
            
            training_X_p1 = X_matrix(sp_tr_1_start: sp_tr_1_end, :);
            training_X_p2 = X_matrix(sp_tr_2_start: sp_tr_2_end, :);
            training_X_p3 = X_matrix(sp_tr_3_start: sp_tr_3_end, :);
            training_X_p4 = X_matrix(sp_tr_4_start: h_X_matrix, :);
            training_X = [training_X_p1; training_X_p2; training_X_p3; training_X_p4];
            
            training_Y_p1 = Y_list(sp_tr_1_start: sp_tr_1_end, :);
            training_Y_p2 = Y_list(sp_tr_2_start: sp_tr_2_end, :);
            training_Y_p3 = Y_list(sp_tr_3_start: sp_tr_3_end, :);
            training_Y_p4 = Y_list(sp_tr_4_start: h_X_matrix, :);
            training_Y = [training_Y_p1; training_Y_p2; training_Y_p3; training_Y_p4];
            
            training_set_X_cells{subjectNum}{musicNum} = training_X;
            training_set_Y_cells{subjectNum}{musicNum} = training_Y;
            test_set_X_cells{subjectNum}{musicNum} = test_X;
            test_set_Y_cells{subjectNum}{musicNum} = test_Y;
        end 
    end
    
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


































































































































































