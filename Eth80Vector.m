% Specify the folder containing quaternion matrices
input_folder = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/1/Resized_RGB_Matrices/Quaternion_Matrices';

% Output file for reduced quaternion matrix
output_file = fullfile(input_folder, 'Final_Reduced_Quaternion_Matrix.mat');

% Number of vectors to retain
k = 10; % Adjust as needed

% Get a list of all quaternion .mat files in the folder
quaternion_files = dir(fullfile(input_folder, '*.mat'));

% Initialize a list to store all quaternion vectors
all_quaternion_vectors = quaternion([],[],[],[]);

% Loop through each quaternion file to form the quaternion matrix
for q = 1:length(quaternion_files)
    % Load the quaternion matrix
    file_path = fullfile(quaternion_files(q).folder, quaternion_files(q).name);
    load(file_path, 'quaternion_matrix');
    real_vector = quaternion_matrix.real(:);
    i_vector = quaternion_matrix.i(:);
    j_vector = quaternion_matrix.j(:);
    k_vector = quaternion_matrix.k(:);
    % Flatten the quaternion matrix into a single vector
quaternion_vector = quaternion(real_vector, i_vector, j_vector, k_vector);

    % Append the quaternion vector as a column
    all_quaternion_vectors = [all_quaternion_vectors quaternion_vector]; % Concatenate column-wise
end

