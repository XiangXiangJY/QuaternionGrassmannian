% Base directory containing subfolders (1, 2, ..., 8)
base_directory = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master';

% List of main directories to process (1 for apples, 2-8 for other categories)
main_dirs = {'1', '2', '3', '4', '5', '6', '7', '8'};

% Output folder for combined quaternion vectors
output_folder = fullfile(base_directory, 'Combined_Quaternion_Vectors');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Process all main directories
for main_dir_idx = 1:length(main_dirs)
    current_dir = fullfile(base_directory, main_dirs{main_dir_idx});
    subdirs = dir(fullfile(current_dir, '*')); % Subfolders within the current directory
    
    % Filter valid subdirectories (exclude '.' and '..')
    subdirs = subdirs([subdirs.isdir] & ~startsWith({subdirs.name}, '.'));
    
    % Process each subdirectory (e.g., Apple_1, Car_1, etc.)
    for sub_idx = 1:length(subdirs)
        subdir_path = fullfile(current_dir, subdirs(sub_idx).name);
        image_extensions = {'*.jpg', '*.png', '*.bmp', '*.tif'};
        image_files = [];
        
        % Collect image files
        for ext = image_extensions
            image_files = [image_files; dir(fullfile(subdir_path, ext{1}))];
        end
        
        if isempty(image_files)
            fprintf('No images found in: %s\n', subdir_path);
            continue;
        end
        
        % Initialize combined quaternion vector
        combined_quaternion_vector = quaternion([], [], [], []);
        
        % Process each image file
        for k = 1:length(image_files)
            image_path = fullfile(image_files(k).folder, image_files(k).name);
            img = imread(image_path);
            img_resized = imresize(img, [20, 20]);
            
            % Convert to quaternion matrix
            quaternion_matrix = quaternion(...
                zeros(size(img_resized, 1), size(img_resized, 2)), ...
                double(img_resized(:, :, 1)), ...
                double(img_resized(:, :, 2)), ...
                double(img_resized(:, :, 3)));
            
            % Flatten into vector and append
            quaternion_vector = quaternion_matrix(:);
            combined_quaternion_vector = [combined_quaternion_vector quaternion_vector];
        end
        
        % Save combined quaternion vector
        category_info = sprintf('%s_%s', main_dirs{main_dir_idx}, subdirs(sub_idx).name);
        output_file = fullfile(output_folder, [category_info, '.mat']);
        save(output_file, 'combined_quaternion_vector', 'category_info');
        fprintf('Processed and saved quaternion vector for %s\n', category_info);
    end
end

disp('All categories processed, quaternion vectors saved.');

% Orthogonalize and compute projector matrices
mat_files = dir(fullfile(output_folder, '*.mat'));
orthogonal_output_directory = fullfile(output_folder, 'Orthogonalized_Matrices');
projector_output_directory = fullfile(output_folder, 'Projector_Matrices');

if ~exist(orthogonal_output_directory, 'dir'), mkdir(orthogonal_output_directory); end
if ~exist(projector_output_directory, 'dir'), mkdir(projector_output_directory); end
t=9;
for k = 1:length(mat_files)
    % Load combined quaternion vector
    mat_file_path = fullfile(mat_files(k).folder, mat_files(k).name);
    data = load(mat_file_path);
    quaternion_matrix0 = data.combined_quaternion_vector;

    [U, S, V] = svd(quaternion_matrix0);

% Retain only the top 9 singular values
Uk = U(:, 1:t);      % n x k matrix
Sk = S(1:t, 1:t);    % k x k diagonal matrix
Vk = V(:, 1:t);  
% Reconstruct the reduced-rank quaternion matrix
quaternion_matrix = mtimes(quaternion_matrix0,Vk);

    [n, m] = size(quaternion_matrix);

    % Initialize orthogonalized matrix
    orthogonal_matrix = quaternion(zeros(n, m), zeros(n, m), zeros(n, m), zeros(n, m));
    q_1 = QuaternionReal(quaternion_matrix(:, 1));
    q_1 = q_1 / normQ(q_1);
    orthogonal_matrix(:, 1) = RealQuaternion(q_1);

    % Gram-Schmidt Process
    for i = 2:m
        q_i = QuaternionReal(quaternion_matrix(:, i));
        for j = 1:i-1
            q_j = QuaternionReal(orthogonal_matrix(:, j));
            projection = timesQ(q_j, timesQ(transQ(q_j), q_i));
            q_i = q_i - projection; % Remove projection
        end
        q_i_norm = normQ(q_i);
        if q_i_norm > 0
            q_i = q_i * (1 / q_i_norm);
        else
            error('Zero norm encountered. Check input vectors for linear dependence.');
        end
        orthogonal_matrix(:, i) = RealQuaternion(q_i);
    end

    % Compute projector matrix
    projector_matrix = timesQ(QuaternionReal(orthogonal_matrix), transQ(QuaternionReal(orthogonal_matrix)));

    % Save orthogonalized and projector matrices
    [~, name, ~] = fileparts(mat_files(k).name);
    save(fullfile(orthogonal_output_directory, [name, '_Orthogonalized.mat']), 'orthogonal_matrix');
    save(fullfile(projector_output_directory, [name, '_Projector.mat']), 'projector_matrix');
    fprintf('Processed and saved matrices for: %s\n', name);
end

disp('Orthogonalization and projector matrices saved.');

% Compute distances between projector matrices
projector_files = dir(fullfile(projector_output_directory, '*_Projector.mat'));
num_files = length(projector_files);
distance_matrix = zeros(num_files);

% Load projector matrices into memory
projector_matrices = cell(num_files, 1);
file_names = cell(num_files, 1);
for k = 1:num_files
    data = load(fullfile(projector_files(k).folder, projector_files(k).name));
    projector_matrices{k} = data.projector_matrix;
    file_names{k} = projector_files(k).name;
end

% Compute distances
[nq, ~] = size(projector_matrices{1});
I = [eye(nq) zeros(nq) zeros(nq) zeros(nq)];
for i = 1:num_files
    for j = i+1:num_files
        P_i = projector_matrices{i};
        P_j = projector_matrices{j};
        diff_matrix = timesQ((I - 2 * P_i), (I - 2 * P_j));
        eigenvalues = quaternion_eigenvalues(RealQuaternion(diff_matrix));
        d = pick_unique_values_with_duplicates(eigenvalues);
        distance_matrix(i, j) = norm(log(d));
        distance_matrix(j, i) = distance_matrix(i, j);
    end
end

% Display distance matrix
disp('Distance Matrix:');
disp(distance_matrix);

% Display pairwise distances
for i = 1:num_files
    for j = i+1:num_files
        fprintf('Distance between %s and %s: %.6f\n', file_names{i}, file_names{j}, distance_matrix(i, j));
    end
end