% ------------ STEP 1: Combine quaternion vectors ------------

% Base folder of image sets
base_directory = '/Users/xiangxiang/Documents/MATLAB/Code fo MTH80 Quaternion Grassmannian Classification/ImageSets';
class_names = {'heavy', 'medium', 'light'};

output_folder = fullfile(base_directory, 'Combined_Quaternion_Vectors');
if ~exist(output_folder, 'dir') 
    mkdir(output_folder);
end

for class_idx = 1:length(class_names) 
    class_dir = fullfile(base_directory, class_names{class_idx});
    video_dirs = dir(class_dir);
    video_dirs = video_dirs([video_dirs.isdir] & ~startsWith({video_dirs.name}, '.'));

    for v = 1:length(video_dirs)
        video_path = fullfile(class_dir, video_dirs(v).name);
        image_files = dir(fullfile(video_path, '*.jpg'));

        if isempty(image_files)
            fprintf('No images in: %s\n', video_path);
            continue;
        end

        combined_quaternion_vector = quaternion([], [], [], []);
        for k = 1:length(image_files)
            img_path = fullfile(image_files(k).folder, image_files(k).name);
            img = imread(img_path);
            img_resized = imresize(img, [24, 24]);

            q_mat = quaternion(...
                zeros(size(img_resized, 1), size(img_resized, 2)), ...
                double(img_resized(:, :, 1)), ...
                double(img_resized(:, :, 2)), ...
                double(img_resized(:, :, 3)));
            q_vec = q_mat(:);
            combined_quaternion_vector = [combined_quaternion_vector q_vec];
        end

        tag = sprintf('%s_%s', class_names{class_idx}, video_dirs(v).name);
        save(fullfile(output_folder, [tag, '.mat']), 'combined_quaternion_vector', 'tag');
        fprintf('Saved quaternion vector for %s\n', tag);
    end
end

disp('All quaternion vectors saved.');

% ------------ STEP 2: PCA + projector matrix ------------

mat_files = dir(fullfile(output_folder, '*.mat'));
orthogonal_output = fullfile(output_folder, 'Orthogonalized_Matrices');
projector_output = fullfile(output_folder, 'Projector_Matrices');
if ~exist(orthogonal_output, 'dir'), mkdir(orthogonal_output); end
if ~exist(projector_output, 'dir'), mkdir(projector_output); end

t = 9;
for k = 1:length(mat_files)
    file_path = fullfile(mat_files(k).folder, mat_files(k).name);
    data = load(file_path);
    quaternion_matrix0 = data.combined_quaternion_vector;

    quaternion_matrix = QuaternionPCA_Column(quaternion_matrix0, t);
    [n, m] = size(quaternion_matrix);

    orthogonal_matrix = quaternion(zeros(n, m), zeros(n, m), zeros(n, m), zeros(n, m));
    q_1 = QuaternionReal(quaternion_matrix(:, 1));
    q_1 = q_1 / normQ(q_1);
    orthogonal_matrix(:, 1) = RealQuaternion(q_1);

    for i = 2:m
        q_i = QuaternionReal(quaternion_matrix(:, i));
        for j = 1:i-1
            q_j = QuaternionReal(orthogonal_matrix(:, j));
            projection = timesQ(q_j, timesQ(transQ(q_j), q_i));
            q_i = q_i - projection;
        end
        q_i_norm = normQ(q_i);
        if q_i_norm > 0
            q_i = q_i * (1 / q_i_norm);
        else
            error('Zero norm encountered.');
        end
        orthogonal_matrix(:, i) = RealQuaternion(q_i);
    end

    projector_matrix = timesQ(QuaternionReal(orthogonal_matrix), transQ(QuaternionReal(orthogonal_matrix)));
    [~, name, ~] = fileparts(mat_files(k).name);
    save(fullfile(orthogonal_output, [name, '_Orthogonalized.mat']), 'orthogonal_matrix');
    save(fullfile(projector_output, [name, '_Projector.mat']), 'projector_matrix');
    fprintf('Projector saved: %s\n', name);
end

% ------------ STEP 3: Compute distances ------------

projector_files = dir(fullfile(projector_output, '*_Projector.mat'));
num_files = length(projector_files);
Data_dir = '/Users/xiangxiang/Documents/MATLAB/Code fo MTH80 Quaternion Grassmannian Classification/hightway';

% Load saved .mat files


load(fullfile(Data_dir, 'distance_matrix.mat'), 'distance_matrix');

%distance_matrix = zeros(num_files);
projector_matrices = cell(num_files, 1);
file_names = cell(num_files, 1);
for k = 1:num_files
    data = load(fullfile(projector_files(k).folder, projector_files(k).name));
    projector_matrices{k} = data.projector_matrix;
    file_names{k} = projector_files(k).name;
end

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
            save(fullfile(Data_dir, 'distance_matrix.mat'), 'distance_matrix');
            fprintf('Checkpoint saved at i = %d, j = %d\n', i, j);
        end
    end
end

save(fullfile(data_dir, 'distance_matrix.mat'), 'distance_matrix');
disp('Distance matrix computation completed and saved.');