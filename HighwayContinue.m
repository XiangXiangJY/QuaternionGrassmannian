% ------------ Resume Distance Matrix Calculation from Separate Files ------------

% Full path to your folder
data_dir = '/Users/xiangxiang/Documents/MATLAB/Code fo MTH80 Quaternion Grassmannian Classification/hightway';

% Load saved .mat files

load(fullfile(data_dir, 'projector_matrix.mat')); 


%load(fullfile(data_dir, 'projector_matrix.mat'), 'projector_matrices');
load(fullfile(data_dir, 'file_names.mat'), 'file_names');
load(fullfile(data_dir, 'distance_matrix.mat'), 'distance_matrix');

num_files = length(projector_);
[nq, ~] = size(projector_matrices{1});
I = [eye(nq) zeros(nq) zeros(nq) zeros(nq)];

% Find the first uncalculated (i, j)
start_i = 1;
start_j = 2;
found = false;

for i = 1:num_files
    for j = i+1:num_files
        if distance_matrix(i, j) == 0
            start_i = i;
            start_j = j;
            found = true;
            break;
        end
    end
    if found
        break;
    end
end

fprintf('Resuming from i = %d, j = %d\n', start_i, start_j);

for i = start_i:num_files
    j_start = (i == start_i) * start_j + (i ~= start_i) * (i + 1);
    for j = j_start:num_files
        if distance_matrix(i, j) ~= 0
            continue;
        end

        progress = ((i - 1) * num_files + j) / (num_files * (num_files - 1) / 2);
        fprintf('Progress: %.2f%% — Computing distance between %d and %d\n', progress * 100, i, j);

        P_i = projector_matrices{i};
        P_j = projector_matrices{j};
        diff_matrix = timesQ((I - 2 * P_i), (I - 2 * P_j));
        eigenvalues = quaternion_eigenvalues(RealQuaternion(diff_matrix));
        d = pick_unique_values_with_duplicates(eigenvalues);
        distance_matrix(i, j) = norm(log(d));
        distance_matrix(j, i) = distance_matrix(i, j);

        if mod(j, 5) == 0
            save(fullfile(data_dir, 'distance_matrix.mat'), 'distance_matrix');
            fprintf('Checkpoint saved at i = %d, j = %d\n', i, j);
        end
    end
end

save(fullfile(data_dir, 'distance_matrix.mat'), 'distance_matrix');
disp('Distance matrix computation completed and saved.');