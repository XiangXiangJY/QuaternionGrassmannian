% Base directory containing subfolders (1, 2, ..., 10)
base_directory = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1';

% List of subdirectories to process
subdirs = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};

% Output folder for combined quaternion vectors with apple information
output_folder = fullfile(base_directory, 'Combined_Quaternion_Vectors_with_Apple');

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
    apple_info = ['Apple_', subdirs{s}]; % E.g., "Apple_1", "Apple_2", ...
    
    % Save the combined quaternion vector and apple information
    output_file = fullfile(output_folder, ['Combined_Quaternion_Vector_Apple_', subdirs{s}, '.mat']);
    save(output_file, 'combined_quaternion_vector', 'apple_info');
    
    fprintf('Processed and saved combined quaternion vector for %s\n', apple_info);
end

disp('All subdirectories have been processed, and quaternion vectors with apple information have been saved.');


%%
% Base directory containing the .mat files
base_directory = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/Combined_Quaternion_Vectors_with_Apple';

% Get a list of all .mat files in the directory
mat_files = dir(fullfile(base_directory, 'Combined_Quaternion_Vector_Apple_*.mat'));

% Check if there are files to process
if isempty(mat_files)
    error('No .mat files found in the specified directory.');
end

% Output directory for the orthogonalized matrices
output_directory = fullfile(base_directory, 'Orthogonalized_Matrices');
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

% Loop through each file and process
for k = 1:length(mat_files)
    % Load the current .mat file
    mat_file_path = fullfile(mat_files(k).folder, mat_files(k).name);
    data = load(mat_file_path);

    % Assume the quaternion matrix is stored in a variable named 'combined_quaternion_vector'
    quaternion_matrix = data.combined_quaternion_vector;

    % Get the number of columns (m) and rows (n) in the matrix
    [n, m] = size(quaternion_matrix);

    
% Initialize the orthogonalized matrix
orthogonal_matrix = quaternion(zeros(n,m),zeros(n,m),zeros(n,m),zeros(n,m));
q_1=QuaternionReal(quaternion_matrix(:, 1));
q_1=q_1/normQ(q_1);
orthogonal_matrix(:,1)=RealQuaternion(q_1);

% Gram-Schmidt Process for Orthogonalization
for i = 2:m
    % Start with the current column vector
    q_i = QuaternionReal(quaternion_matrix(:, i));
    
    
    % Subtract projections onto already orthogonalized vectors
    for j = 1:i-1
        q_j = QuaternionReal(orthogonal_matrix(:, j));
projection =timesQ(q_j,timesQ(transQ(q_j),q_i)); % Projection of q_i onto q_j
        q_i = q_i - projection; % Remove the projection component
    end
    
    % Normalize the orthogonalized vector
    q_i_norm = normQ(q_i);
     if q_i_norm > 0
        q_i = q_i * (1 / q_i_norm);
    else
        error('Zero norm encountered. Check input vectors for linear dependence.');
    end
   
    
   
    orthogonal = true; % Flag to track orthogonality
    for j = 1:i-1
        q_j = QuaternionReal(orthogonal_matrix(:, j));
        dot_product = normQ((timesQ(transQ(q_j), q_i))); % Real part of the quaternion dot product
        
        if abs(dot_product) > 1e-6 % Tolerance for numerical precision
            orthogonal = false;
            fprintf('Warning: Vector %d is not orthogonal to Vector %d (Dot Product: %e)\n', i, j, dot_product);
        end
    end
    

    % Store the orthogonalized vector in the output matrix
    if orthogonal
        orthogonal_matrix(:, i) = RealQuaternion(q_i);
    else
        fprintf('Skipping Vector %d due to orthogonality issues.\n', i);
    end
end

    %Take ortgogonal_matrix*ortgogonal_matrix^*, get one point (projector)
    %in grassmannian
       %MO=QuaternionReal(orthogonal_matrix);
       %P=timesQ(MO,transQ(MO));
       %Point_Grass=RealQuaternion(P);

    % Save the orthogonalized matrix
    [~, name, ~] = fileparts(mat_files(k).name); % Extract file name without extension
    save_path = fullfile(output_directory, [name, '_Orthogonalized.mat']);
    save(save_path, 'orthogonal_matrix');
    %save_path = fullfile(output_directory, [name, '_Point.mat']);
    %save(save_path, 'Point_Grass');


    fprintf('Processed and saved orthogonalized matrix for file: %s\n', mat_files(k).name);
end

disp('All files have been processed and orthogonalized matrices saved.');

%%
% get grassmannian points corresponding to image sets

% Base directory containing the .mat files
base_directory = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/Combined_Quaternion_Vectors_with_Apple';

% Get a list of all .mat files in the directory
mat_files = dir(fullfile(base_directory, 'Combined_Quaternion_Vector_Apple_*.mat'));

% Check if there are files to process
if isempty(mat_files)
    error('No .mat files found in the specified directory.');
end

% Output directories for the orthogonalized matrices and projector matrices
orthogonal_output_directory = fullfile(base_directory, 'Orthogonalized_Matrices');
projector_output_directory = fullfile(base_directory, 'Projector_Matrices');

% Create the directories if they don't exist
if ~exist(orthogonal_output_directory, 'dir')
    mkdir(orthogonal_output_directory);
end
if ~exist(projector_output_directory, 'dir')
    mkdir(projector_output_directory);
end

% Loop through each file and process
for k = 1:length(mat_files)
    % Load the current .mat file
    mat_file_path = fullfile(mat_files(k).folder, mat_files(k).name);
    data = load(mat_file_path);

    % Assume the quaternion matrix is stored in a variable named 'combined_quaternion_vector'
    quaternion_matrix = data.combined_quaternion_vector;

    % Get the number of columns (m) and rows (n) in the matrix
    [n, m] = size(quaternion_matrix);

    % Initialize the orthogonalized matrix
    orthogonal_matrix = quaternion(zeros(n, m), zeros(n, m), zeros(n, m), zeros(n, m));
    q_1 = QuaternionReal(quaternion_matrix(:, 1));
    q_1 = q_1 / normQ(q_1);
    orthogonal_matrix(:, 1) = RealQuaternion(q_1);

    % Gram-Schmidt Process for Orthogonalization
    for i = 2:m
        % Start with the current column vector
        q_i = QuaternionReal(quaternion_matrix(:, i));
        
        % Subtract projections onto already orthogonalized vectors
        for j = 1:i-1
            q_j = QuaternionReal(orthogonal_matrix(:, j));
            projection = timesQ(q_j, timesQ(transQ(q_j), q_i)); % Projection of q_i onto q_j
            q_i = q_i - projection; % Remove the projection component
        end
        
        % Normalize the orthogonalized vector
        q_i_norm = normQ(q_i);
        if q_i_norm > 0
            q_i = q_i * (1 / q_i_norm);
        else
            error('Zero norm encountered. Check input vectors for linear dependence.');
        end

        % Check orthogonality with all previous vectors
        orthogonal = true; % Flag to track orthogonality
        for j = 1:i-1
            q_j = QuaternionReal(orthogonal_matrix(:, j));
            dot_product = normQ(timesQ(transQ(q_j), q_i)); % Real part of the quaternion dot product
            
            if abs(dot_product) > 1e-6 % Tolerance for numerical precision
                orthogonal = false;
                fprintf('Warning: Vector %d is not orthogonal to Vector %d (Dot Product: %e)\n', i, j, dot_product);
            end
        end

        % Store the orthogonalized vector in the output matrix
        if orthogonal
            orthogonal_matrix(:, i) = RealQuaternion(q_i);
        else
            fprintf('Skipping Vector %d due to orthogonality issues.\n', i);
        end
    end

    % Compute the projector quaternion matrix
    apple41 = timesQ(QuaternionReal(orthogonal_matrix), transQ(QuaternionReal(orthogonal_matrix)));

    % Save the orthogonalized matrix in the orthogonal directory
    [~, name, ~] = fileparts(mat_files(k).name); % Extract file name without extension
    save_path_orthogonal = fullfile(orthogonal_output_directory, [name, '_Orthogonalized.mat']);
    save(save_path_orthogonal, 'orthogonal_matrix');

    % Save the projector quaternion matrix in the projector directory
    save_path_projector = fullfile(projector_output_directory, [name, '_Projector.mat']);
    save(save_path_projector, 'apple41');

    fprintf('Processed and saved orthogonalized and projector matrices for file: %s\n', mat_files(k).name);
end

disp('All files have been processed, and matrices saved.');

% Directory containing the projector matrices
projector_directory = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/Combined_Quaternion_Vectors_with_Apple/Projector_Matrices';

% Get a list of all .mat files in the projector directory
projector_files = dir(fullfile(projector_directory, '*_Projector.mat'));

% Check if there are files to process
if isempty(projector_files)
    error('No projector files found in the specified directory.');
end

% Load all projector matrices into a cell array
num_files = length(projector_files);
projector_matrices = cell(num_files, 1);
file_names = cell(num_files, 1);

for k = 1:num_files
    % Load the current projector file
    projector_file_path = fullfile(projector_files(k).folder, projector_files(k).name);
    data = load(projector_file_path);
    
    % Assume the projector quaternion matrix is stored in a variable named 'apple41'
    projector_matrices{k} = data.apple41;
    file_names{k} = projector_files(k).name;
end

% Initialize a matrix to store distances
distance_matrix = zeros(num_files, num_files);
[nq,mq]=size(projector_matrices{1});
I=[eye(nq) zeros(nq) zeros(nq) zeros(nq)];

% Compute distances between all pairs of projector matrices
for i = 1:num_files
    for j = i+1:num_files
        % Compute the Frobenius norm-based distance between two quaternion matrices
        diff_matrix = timesQ((I-2.*projector_matrices{i}), (I-2.*projector_matrices{j}));
        D=quaternion_eigenvalues(RealQuaternion(diff_matrix));
        d=pick_unique_values_with_duplicates(D);
        distance_matrix(i, j) =norm(log(d)); % Frobenius norm of the difference
        distance_matrix(j, i) = distance_matrix(i, j); % Symmetric
    end
end

% Display the distance matrix
disp('Distance Matrix:');
disp(distance_matrix);

% Display pairwise distances with file names
for i = 1:num_files
    for j = i+1:num_files
        fprintf('Distance between %s and %s: %.6f\n', file_names{i}, file_names{j}, distance_matrix(i, j));
    end
end