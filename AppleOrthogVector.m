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