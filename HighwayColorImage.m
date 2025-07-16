% Paths
convertedFolder = '/Users/xiangxiang/Desktop/Revise JMIV/video/converted';
outputRoot = fullfile(pwd, 'ImageSets');

% Load labels
load('imagemaster.mat');  % expects imagemaster{i}.root and .class

% Create class folders
classes = {'heavy', 'medium', 'light'};
for i = 1:numel(classes)
    classPath = fullfile(outputRoot, classes{i});
    if ~exist(classPath, 'dir')
        mkdir(classPath);
    end
end

% Process each video
for i = 1:length(imagemaster)
    videoName = imagemaster{i}.root;
    className = imagemaster{i}.class;

    inputVideo = fullfile(convertedFolder, [videoName, '_mjpeg.avi']);
    if ~isfile(inputVideo)
        fprintf(' Skipping (not found): %s\n', inputVideo);
        continue;
    end

    % Output folder: ImageSets/class/videoName/
    outputFolder = fullfile(outputRoot, className, videoName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    try
        v = VideoReader(inputVideo);
    catch
        fprintf(' Error opening video: %s\n', inputVideo);
        continue;
    end

    frameIdx = 1;
    while hasFrame(v)
        frame = readFrame(v);  % RGB image
        filename = fullfile(outputFolder, sprintf('frame_%04d.jpg', frameIdx));
        imwrite(frame, filename);
        frameIdx = frameIdx + 1;
    end

    fprintf('Saved %d frames for %s (%s)\n', frameIdx-1, videoName, className);
end

fprintf('\n All video frames saved into class-based image sets under:\n%s\n', outputRoot);