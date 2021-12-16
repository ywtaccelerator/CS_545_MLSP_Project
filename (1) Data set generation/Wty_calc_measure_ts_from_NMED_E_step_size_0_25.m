% This script has not  been tested on cloud machine
% Workspace directory: 
% /home/soren/Documents/Wentao/code/GCP_matlab_workspace/11. ML for NMED-E/20210823_Calculating neural measure time series from the EEG data


% ======== Please specify the following information for calculation ========
% Please specify the window size in second:
windowSize = 3;

% Please specify the overlaping size (in seconds) of sliding window:
% overlapSize = 2; % use an overlap of 2/3 window size
overlapSize = 2.75; % currently use an overlap of 2.75 seconds,
                   % i.e. the step size is 0.25 second.
                   % Thus the sampling rate of result measure time series
                   % is 1/0.25 = 4 Hz

% Please specify the sampling rate of MEG in Hz:
samplingRate = 125;

% Output directory
outputFolderName = 'NMED_E_measure_time_series_output'; % output folder name
mkdir(outputFolderName); 

% Input EEG time series are stored in 'my_EEG_files'

path_data_file_122 = ['/home/soren/Documents/Wentao/data/NMED-E/EEG/', ...
             'CleanEEG_stim122.mat'];
path_data_file_22 = ['/home/soren/Documents/Wentao/data/NMED-E/EEG/', ...
             'CleanEEG_stim22.mat'];
path_data_file_123 = ['/home/soren/Documents/Wentao/data/NMED-E/EEG/', ...
             'CleanEEG_stim123.mat'];      
path_data_file_23 = ['/home/soren/Documents/Wentao/data/NMED-E/EEG/', ...
             'CleanEEG_stim23.mat'];


% ======== Please specify the above information for calculation ============

load(path_data_file_122);
load(path_data_file_22);
load(path_data_file_123);
load(path_data_file_23);

eeg_data_set{1} = eeg122;
eeg_data_set{2} = eeg22;
eeg_data_set{3} = eeg123;
eeg_data_set{4} = eeg23;

clear eeg122 eeg22 eeg123 eeg23;

% files = dir(fullfile(pathInput));


% % convert the window size into number of samples
% windowLength = windowSize * samplingRate;
% 
% % calculate the step length (in number of samples)
% stepLength = round( (windowSize - overlapSize) * samplingRate );

stepSize = windowSize - overlapSize;

measureTimeSeriesFromEEG = {};
% subjectNumber = 1;

% load('../MEG_intersection_list_of_subject_numbers.mat');


for subjectNum = 1:23 % 23 subjects in total
%     fprintf('Subject # %d\n\n', subjectNum);

    measureTimeSeriesOfTheTrial = {};
    for dataSetNum = 1:4
        % for each subject
        eeg_data = eeg_data_set{dataSetNum};
        
        table0 = eeg_data(:,:,subjectNum);
        [numChannels, numSamples] = size(table0);
        
        % creating a timeLine
        timeLine = [];
        timeStamp = 0; % a time stamp starts from 0
        stepLength = 1/samplingRate;
        for sampleNum = 1:numSamples
            timeLine = [timeLine, timeStamp];
            timeStamp = timeStamp + stepLength;
        end
        
        
        measureTimeSeriesOfTheChannel = {};

        % calculate measures for each channel's time series
        parfor channelNum = 1:numChannels
            timeSeries = table0(channelNum,:)';

            % temporarily replace all NaNs by 0s
            timeSeries(find(isnan(timeSeries)==1)) = 0;

            resultMF = [];
            resultSpectralEntropy = [];
            resultDeltaBandPower = [];
            resultThetaBandPower = [];
            resultGammaBandPower = [];
            resultBetaBandPower = [];
            resultAlphaBandPower = [];
            resultACW = [];
            resultPLE = [];
            resultSD = [];
            resultLZC = [];
            resultShannonEntropy = [];

            endTimeStampOfTheSubject = timeLine(length(timeLine));

            currentWindowNum = 1; % start from the 1st window
            % for each of the (numWindows - 1) whole windows
            startTimeStamp = (currentWindowNum - 1) * stepSize;
            endTimeStamp = startTimeStamp + windowSize;

            while endTimeStamp + stepSize <= endTimeStampOfTheSubject
                for i = 1:length(timeLine)
                    if timeLine(i) >= startTimeStamp
                        startSampleNumber = i;
                        break;
                    end
                end
                for j = startSampleNumber:length(timeLine)
                    if timeLine(j) > endTimeStamp
                        endSampleNumber = j - 1;
                        break;
                    end
                end

    %                 startSampleNumber = 1 + (windowNumber - 1) * stepLength;
    %                 endSampleNumber = startSampleNumber + windowLength - 1;
                windowSamples = timeSeries(startSampleNumber : endSampleNumber);

                MF = medfreq(windowSamples, samplingRate,[1 40]);
                resultMF = [resultMF; MF];

                % !!!!! currently cannot be calculated using the MATLAB installed on
                % the Cruncher,
                % hence, currently let us fill the column for spectralEntropy
                % by all 0's
    %             spectralEntropy = pentropy(windowSamples, samplingRate, 'Instantaneous', false);
    %             resultSpectralEntropy = [resultSpectralEntropy; spectralEntropy];
                resultSpectralEntropy = [resultSpectralEntropy; 0];


                deltaBandPower = bandpower(windowSamples, samplingRate, [1, 4]);
                resultDeltaBandPower = [resultDeltaBandPower; deltaBandPower];

                thetaBandPower = bandpower(windowSamples, samplingRate, [4, 8]);
                resultThetaBandPower = [resultThetaBandPower; thetaBandPower];

                % changed the bandpass for calculating gamma band power,
                % from [30, 50] to [30, 40], to be consistent with the
                % calculation of neural measure time series on the music side.
                % 
                % Because when calculating Gamma Band Power on the music side,
                % due to the sampling rate of short-term acoustic feature time series is
                % 1s / (25ms/2) = 80 Hz, thus the upper limit for calculating
                % bandpower is 40. So if we run the following code:
                % gammaBandPower = bandpower(windowSamples, samplingRate, [30, 50]);
                % there will be an error saying:
                % "The frequency range must be within the range of F for the specified
                % input."
                % Therefore, I restricted the upper limit of bandpass to
                % 40, i.e., calculate the Gamma Band Power using bandpass [30, 40],
                % instead of using [30, 50].
                gammaBandPower = bandpower(windowSamples, samplingRate, [30, 40]);
                resultGammaBandPower = [resultGammaBandPower; gammaBandPower];

                betaBandPower = bandpower(windowSamples, samplingRate, [13, 30]);
                resultBetaBandPower = [resultBetaBandPower; betaBandPower];

                alphaBandPower = bandpower(windowSamples, samplingRate, [8, 13]);
                resultAlphaBandPower = [resultAlphaBandPower; alphaBandPower];


                % calculating ACW
                % There are some problems with the use of ACW_estimation on
                % window samples, I cannot calculate ACW on the whole window of sample
                %%%%%% The following code is modified from ACW_estimation.m %%%%%%
                EEG = windowSamples;
                srate = samplingRate;

                overlap=50;
                lag=(length(EEG)-1)/srate;
                lag=floor(lag*srate);

                % ACF computation into the window (normalized between -1 and 1)
    %                 ACF(1,:)=xcorr(EEG,lag,'coeff'); % changed the way of
    %                 calculating ACF to the following way to enable parallel
    %                 processing
                ACF = xcorr(EEG,lag,'coeff')';
                % As in Honey et al. 2012, ACF is averaged
                ACF_mean=mean(ACF,1);
                % Look for the index where the mean of the ACF is maximum
                [~,index]=max(ACF_mean);
                % Number of indices over the half in the right side of the lobe
                my_index=ACF_mean>=max(ACF_mean)/2;
                myindex=my_index(index:end);
                under_half=find(myindex==0);
                first_under_half=under_half(1);

                ACW_samples=2*(first_under_half-1)-1;   % ACW in samples
                ACW=ACW_samples/srate;                     % ACW in time
                %%%%%% End of the modification from ACW_estimation.m %%%%%%
                resultACW = [resultACW; ACW];


                PLE = JF_power_law(windowSamples, 1/samplingRate, 1, 45); %45 good for EEG,but not very much for the music side
                resultPLE = [resultPLE; PLE];

                SD = std(windowSamples);
                resultSD = [resultSD; SD];

                LZC = LZC_estimation(windowSamples);
                resultLZC = [resultLZC; LZC];

                shanEn = wentropy(windowSamples,'shannon')/length(windowSamples);
                resultShannonEntropy = [resultShannonEntropy; shanEn];


                currentWindowNum = currentWindowNum + 1;
                startTimeStamp = (currentWindowNum - 1) * stepSize;
                endTimeStamp = startTimeStamp + windowSize;
            end

            columns = { 'resultMF', ...
                        'resultSpectralEntropy', ...
                        'resultDeltaBandPower', ...
                        'resultThetaBandPower', ...
                        'resultGammaBandPower', ...
                        'resultBetaBandPower', ...
                        'resultAlphaBandPower', ...
                        'resultACW', ...
                        'resultPLE', ...
                        'resultSD', ...
                        'resultLZC', ...
                        'resultShannonEntropy'};

    %         columns = { 'resultMF', ...
    %                     'resultACW', ...
    %                     'resultPLE', ...
    %                     'resultLZC'};

            measureTimeSeriesOfTheChannel{channelNum} = table(  resultMF,...
                                                                resultSpectralEntropy,...
                                                                resultDeltaBandPower,...
                                                                resultThetaBandPower,...
                                                                resultGammaBandPower, ...
                                                                resultBetaBandPower,...
                                                                resultAlphaBandPower,...
                                                                resultACW,...
                                                                resultPLE,...
                                                                resultSD,...
                                                                resultLZC,...
                                                                resultShannonEntropy,...
                                                                'VariableNames', ...
                                                                columns);

    %         measureTimeSeriesOfTheChannel{channelNum} = table(  resultMF,...
    %                                                             resultACW,...
    %                                                             resultPLE,...
    %                                                             resultLZC,...
    %                                                             'VariableNames', ...
    %                                                             columns);                                                

            fprintf('Calculated all measures for:\n - subject #%d - trial #%d - Channel #%d\n\n', subjectNum, dataSetNum, channelNum);             
        end
        
        measureTimeSeriesOfTheTrial{dataSetNum} = measureTimeSeriesOfTheChannel;
    end
    measureTimeSeriesFromEEG{subjectNum} = measureTimeSeriesOfTheTrial;

    %     subjectNumber = subjectNumber + 1;
end

outputPath = [outputFolderName, '/measure_ts_from_NMED_E_step_size_0_25.mat'];
save(outputPath,'measureTimeSeriesFromEEG', '-v7.3');



