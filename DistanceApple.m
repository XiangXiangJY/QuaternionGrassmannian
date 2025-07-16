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
disp('Distance Matrix (Frobenius Norm):');
disp(distance_matrix);

% Display pairwise distances with file names
for i = 1:num_files
    for j = i+1:num_files
        fprintf('Distance between %s and %s: %.6f\n', file_names{i}, file_names{j}, distance_matrix(i, j));
    end
end