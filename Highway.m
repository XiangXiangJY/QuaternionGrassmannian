% 配置路径
originalFolder = '/Users/xiangxiang/Desktop/Revise JMIV/video';
convertedFolder = fullfile(originalFolder, 'converted');
ffmpegPath = '/Users/xiangxiang/Downloads/ffmpeg';  % 你下载的ffmpeg可执行文件路径

% 创建输出文件夹（如不存在）
if ~exist(convertedFolder, 'dir')
    mkdir(convertedFolder);
end

% 加载标签
load('imagemaster.mat');  % 包含 imagemaster{i}.root

% 遍历每个视频
for i = 1:length(imagemaster)
    name = imagemaster{i}.root;

    inputVideo = fullfile(originalFolder, [name, '.avi']);
    outputVideo = fullfile(convertedFolder, [name, '_mjpeg.avi']);

    % 跳过不存在的原视频
    if ~isfile(inputVideo)
        fprintf('❌ Missing file: %s\n', inputVideo);
        continue;
    end

    % 跳过已转码的视频
    if isfile(outputVideo)
        fprintf('✅ Already converted: %s\n', name);
        continue;
    end

    % 构造 FFmpeg 命令
    ffmpegCmd = sprintf('"%s" -y -i "%s" -c:v mjpeg -q:v 3 -an "%s"', ...
        ffmpegPath, inputVideo, outputVideo);

    % 执行
    [status, result] = system(ffmpegCmd);

    % 输出信息
    if status == 0
        fprintf('Converted: %s\n', name);
    else
        fprintf('Failed to convert: %s\nError: %s\n', name, result);
    end
end

fprintf('\n All videos processed. Transcoded files saved in ：\n%s\n', convertedFolder);