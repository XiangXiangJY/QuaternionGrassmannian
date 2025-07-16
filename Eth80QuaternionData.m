% Specify the folder containing resized RGB matrices
input_folder = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/1/Resized_RGB_Matrices';

% Output folder for quaternion matrices
output_folder = fullfile(input_folder, 'Quaternion_Matrices');

% Create the output folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get a list of all .mat files in the folder
mat_files = dir(fullfile(input_folder, '*.mat'));

% Loop through each file
for k = 1:length(mat_files)
    % Load the resized RGB matrix
    file_path = fullfile(mat_files(k).folder, mat_files(k).name);
    load(file_path, 'img_resized');
    
    % Extract RGB channels
    red_channel = double(img_resized(:, :, 1)); % i part
    green_channel = double(img_resized(:, :, 2)); % j part
    blue_channel = double(img_resized(:, :, 3)); % k part
    
    % Form a pure quaternion matrix (real part is zero)
    quaternion_matrix = struct('real', zeros(size(red_channel)), ...
                                'i', red_channel, ...
                                'j', green_channel, ...
                                'k', blue_channel);
    
    % Save the quaternion matrix
    [~, name, ~] = fileparts(mat_files(k).name); % Get the file name without extension
    save_path = fullfile(output_folder, [name, '_Quaternion.mat']);
    save(save_path, 'quaternion_matrix');
    
    fprintf('Processed and saved quaternion matrix for: %s\n', mat_files(k).name);
end

disp('All images have been converted to quaternion matrices and saved.');