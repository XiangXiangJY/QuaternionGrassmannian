% Assume you have a distance matrix named distance_matrix (80x80)
% and it's ordered such that each group of 10 belongs to one class

% Set random seed for reproducibility
a= load('distance_matrixPCA.mat');
distance_matrix = a.distance_matrix;
% Assume distance_matrix is already loaded (80x80)
rng('shuffle');  % For reproducibility

num_classes = 8;
samples_per_class = 10;
selected_per_class = 5;

% --- Step 1: Select 5 random samples per class for training ---
selected_indices = [];
for class_idx = 0:(num_classes - 1)
    class_start = class_idx * samples_per_class + 1;
    class_range = class_start:(class_start + samples_per_class - 1);
    selected = randsample(class_range, selected_per_class);
    selected_indices = [selected_indices, selected];
end

% Sort and label training data
[selected_indices_sorted, sort_order] = sort(selected_indices);
labels_train = repelem(1:num_classes, selected_per_class);
labels_train = labels_train(sort_order);

% Get training submatrix (40x40)
D_train = distance_matrix(selected_indices_sorted, selected_indices_sorted);

% --- Step 2: Pick one random test sample from remaining 40 ---
all_indices = 1:80;
remaining_indices = setdiff(all_indices, selected_indices_sorted);
test_index = randsample(remaining_indices, 1);

% Get distances from test sample to training samples (1x40)
D_test_row = distance_matrix(test_index, selected_indices_sorted);

% Combine into extended 41x41 distance matrix
D_extended = [D_train, D_test_row'; D_test_row, 0];

% --- Step 3: Perform MDS on the 41x41 matrix ---
[Y_ext, stress] = mdscale(D_extended, 2, 'Criterion', 'stress');

Y_train = Y_ext(1:40, :);
Y_test = Y_ext(41, :);

% --- Step 4: Plot ---
% Define colors and markers
colors = lines(num_classes);
markers = {'o', '+', '*', 'x', 's', 'd', '^', 'v'};

figure;

hold on;
for i = 1:num_classes
    idx = labels_train == i;
    scatter(Y_train(idx,1), Y_train(idx,2), 80, ...
        'Marker', markers{i}, ...
        'MarkerEdgeColor', colors(i,:), ...
        'LineWidth', 1.2, ...
        'DisplayName', ['Class ' num2str(i)]);
end

% Plot the test point in red
scatter(Y_test(1), Y_test(2), 100, 'ro', 'filled', ...
    'DisplayName', 'Test Point');

 
%title('MDS Plot with Test Sample from Quaternionic Grassmannian Distance Matrix', 'FontSize', 14);
xlabel('MDS Dimension 1');
ylabel('MDS Dimension 2');
legend('Location', 'bestoutside');
axis equal;