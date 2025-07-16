% Base directory containing subfolders (1, 2, ..., 10)
base_directory = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/2';

% List of subdirectories to process
subdirs = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};

% Output folder for combined quaternion vectors with apple information
output_folder = fullfile(base_directory, 'Combined_Quaternion_Vectors_with_Car');

% Create the output folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Loop through each subdirectory
for s = 1:length(subdirs)
    % Specify the current subdirectory
    subdir_path = fullfile(base_directory, subdirs{s});
    
    % Get a list of all image files in the subdirectory
    image_extensions = {'*.jpg', '*.png', '*.bmp', '*.tif'};
    image_files = [];
    for ext = image_extensions
        image_files = [image_files; dir(fullfile(subdir_path, ext{1}))];
    end
    
    % Check if there are any images to process
    if isempty(image_files)
        fprintf('No images found in subdirectory: %s\n', subdirs{s});
        continue;
    end
    
    % Initialize the combined quaternion vector for this apple
    combined_quaternion_vector = quaternion([], [], [], []);
    
    % Loop through each image file
    for k = 1:length(image_files)
        % Read and resize the image to 20x20
        image_path = fullfile(image_files(k).folder, image_files(k).name);
        img = imread(image_path);
        img_resized = imresize(img, [20, 20]);
        
        % Convert the resized image into a quaternion matrix
        quaternion_matrix = quaternion(...
            zeros(size(img_resized, 1), size(img_resized, 2)), ... % Real part is 0
            double(img_resized(:, :, 1)), ... % i part (red channel)
            double(img_resized(:, :, 2)), ... % j part (green channel)
            double(img_resized(:, :, 3)));   % k part (blue channel)
        
        % Flatten the quaternion matrix into a vector
        quaternion_vector = quaternion_matrix(:);
        
        % Append the quaternion vector to the combined quaternion vector
        combined_quaternion_vector = [combined_quaternion_vector quaternion_vector];
    end
    
    % Add apple information
    car_info = ['Car_', subdirs{s}]; % E.g., "Car_1", "Car_2", ... 
    
    % Save the combined quaternion vector and apple information
    output_file = fullfile(output_folder, ['Combined_Quaternion_Vector_Car_', subdirs{s}, '.mat']);
    save(output_file, 'combined_quaternion_vector', 'car_info');
    
    fprintf('Processed and saved combined quaternion vector for %s\n', car_info);
end

disp('All subdirectories have been processed, and quaternion vectors with car information have been saved.');