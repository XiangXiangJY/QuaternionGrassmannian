% --------- STEP 3 Parallelized with parfor ----------
projector_files = dir(fullfile(projector_output, '*_Projector.mat'));
num_files = length(projector_files);
Data_dir = '/Users/xiangxiang/Documents/MATLAB/Code fo MTH80 Quaternion Grassmannian Classification/hightway';

load(fullfile(Data_dir, 'distance_matrix.mat'), 'distance_matrix');

projector_matrices = cell(num_files, 1);
file_names = cell(num_files, 1);
for k = 1:num_files
    data = load(fullfile(projector_files(k).folder, projector_files(k).name));
    projector_matrices{k} = data.projector_matrix;
    file_names{k} = projector_files(k).name;
end

[nq, ~] = size(projector_matrices{1});
I = [eye(nq) zeros(nq) zeros(nq) zeros(nq)];


% Only compute for i = 59 to 69 and j > i
target_i = 59:69;
remaining_pairs = [];
for i = target_i
    for j = i+1:num_files
        if distance_matrix(i, j) == 0
            remaining_pairs = [remaining_pairs; i, j];
        end
    end
end

fprintf('Selected pairs for i = 59 to 69: %d pairs\n', size(remaining_pairs, 1));

% Start parpool(2) if not already running
if isempty(gcp('nocreate'))
    parpool(2);
end

% Preallocate results
temp_results(size(remaining_pairs, 1)) = struct('i', 0, 'j', 0, 'value', 0);

% Compute distances in parallel
parfor idx = 1:size(remaining_pairs, 1)
    i = remaining_pairs(idx, 1);
    j = remaining_pairs(idx, 2);

    P_i = projector_matrices{i};
    P_j = projector_matrices{j};

    diff_matrix = timesQ((I - 2 * P_i), (I - 2 * P_j));
    eigenvalues = quaternion_eigenvalues(RealQuaternion(diff_matrix));
    d = pick_unique_values_with_duplicates(eigenvalues);
    result = norm(log(d));

    temp_results(idx).i = i;
    temp_results(idx).j = j;
    temp_results(idx).value = result;
end

% Fill values into distance_matrix
for idx = 1:length(temp_results)
    i = temp_results(idx).i;
    j = temp_results(idx).j;
    distance_matrix(i, j) = temp_results(idx).value;
    distance_matrix(j, i) = temp_results(idx).value;
end

% Save updated distance matrix
save(fullfile(data_dir, 'distance_matrix.mat'), 'distance_matrix');
disp('Partial distance matrix updated and saved (i = 59 to 69).');