% Define original and converted video folder paths
originalFolder = '/Users/xiangxiang/Desktop/Revise JMIV/video';
convertedFolder = fullfile(originalFolder, 'converted');
ffmpegPath = '/Users/xiangxiang/Downloads/ffmpeg';  % Path to ffmpeg executable

% Create output folder if it does not exist
if ~exist(convertedFolder, 'dir')
    mkdir(convertedFolder);
end

% Load label data
load('imagemaster.mat');  % expects imagemaster{i}.root

% Process each video
for i = 1:length(imagemaster)
    name = imagemaster{i}.root;

    inputVideo = fullfile(originalFolder, [name, '.avi']);
    outputVideo = fullfile(convertedFolder, [name, '_mjpeg.avi']);

    % Skip missing files
    if ~isfile(inputVideo)
        fprintf('Missing file: %s\n', inputVideo);
        continue;
    end

    % Skip if already converted
    if isfile(outputVideo)
        fprintf('Already converted: %s\n', name);
        continue;
    end

    % Build FFmpeg command
    ffmpegCmd = sprintf('"%s" -y -i "%s" -c:v mjpeg -q:v 3 -an "%s"', ...
        ffmpegPath, inputVideo, outputVideo);

    % Run the command
    [status, result] = system(ffmpegCmd);

    % Print result
    if status == 0
        fprintf('Converted: %s\n', name);
    else
        fprintf(' Failed to convert: %s\nError: %s\n', name, result);
    end
end

fprintf('\n All videos have been processed.\nConverted videos are saved in:\n%s\n', convertedFolder);