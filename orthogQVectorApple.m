% Specify the path to your .mat file
mat_file = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/Combined_Quaternion_Vectors_with_Apple/Combined_Quaternion_Vector_Apple_1.mat';

%
% Specify the path to your .mat file
%mat_file = '/path/to/your/matfile.mat';

% Load the data from the .mat file
data = load(mat_file);

% Assume the quaternion matrix is stored in a variable named 'quaternion_matrix'


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

% Save the orthogonalized matrix
save_path = '/Users/xiangxiang/Library/CloudStorage/OneDrive-UniversityofNevada,Reno/2024/Fall 2024/Math 799/Quternion Grassmannian in Image set reconigtion/Dataset/ETH-80-master/ETH-80-master/1/Combined_Quaternion_Vectors_with_Apple/orthogonalized_matrix1.mat';
save(save_path, 'orthogonal_matrix');

disp(['Orthogonalized matrix saved to: ', save_path]);

