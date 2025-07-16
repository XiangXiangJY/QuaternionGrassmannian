% Specify the path to your OneDrive folder
onedrive_folder = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/1';

% Output folder to save resized RGB matrices
output_folder = fullfile(onedrive_folder, 'Resized_RGB_Matrices'); 

% Create the output folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get a list of all image files in the folder (common image extensions)
image_extensions = {'*.jpg', '*.png', '*.bmp', '*.tif'}; 
image_files = [];
for ext = image_extensions
    image_files = [image_files; dir(fullfile(onedrive_folder, ext{1}))];
end

% Check if there are any images to process
if isempty(image_files)
    disp('No images found in the specified folder.');
    return;
end

% Desired size for resizing
resize_size = [20, 20]; % 20x20 pixels

% Loop through each image file
for k = 1:length(image_files)
    % Get the full file path
    image_path = fullfile(image_files(k).folder, image_files(k).name);
    
    % Read the image
    img = imread(image_path);
    
    % Resize the image to 20x20
    img_resized = imresize(img, resize_size);
    
    % Get the file name without extension
    [~, name, ~] = fileparts(image_files(k).name);
    
    % Save the resized RGB matrix
    save_path = fullfile(output_folder, [name, '_Resized_RGB.mat']);
    save(save_path, 'img_resized');
    
    fprintf('Processed, resized, and saved RGB matrix for: %s\n', image_files(k).name);
end

disp('All images have been resized and RGB matrices saved.');