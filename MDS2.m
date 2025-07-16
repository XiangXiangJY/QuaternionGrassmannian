% Load distance matrix
a = load('highwaydistance_matrix.mat');
distance_matrix = a.distance_matrix;

% Assume distance_matrix is 254x254, already in workspace
% Class order: 1–44, 45–209, 210–254
rng('shuffle');  % Use 'rng(1)' for fixed results

total_samples = 254;
train_size = 195;

all_indices = 1:total_samples;

% Step 1: Randomly select 195 samples for training
train_indices = randsample(all_indices, train_size);

% Step 2: Select 1 test sample from the remaining 59
remaining_indices = setdiff(all_indices, train_indices);
test_index = randsample(remaining_indices, 1);

% Step 3: Build full 196x196 distance matrix for MDS
selected_indices = [train_indices, test_index];
D_sub = distance_matrix(selected_indices, selected_indices);

% Step 4: Perform MDS
[Y, ~] = mdscale(D_sub, 2, 'Criterion', 'stress');

Y_train = Y(1:end-1, :);
Y_test = Y(end, :);

% Step 5: Assign class labels for coloring
labels = zeros(train_size,1);
for i = 1:train_size
    idx = train_indices(i);
    if idx <= 44
        labels(i) = 1;
    elseif idx <= 209
        labels(i) = 2;
    else
        labels(i) = 3;
    end
end

% Step 6: Simple plot
figure;
gscatter(Y_train(:,1), Y_train(:,2), labels, ['b','g','m'], 'o', 8);
hold on;
scatter(Y_test(1), Y_test(2), 100, 'r', 'filled', 'DisplayName', 'Test Point');
%title('MDS Plot of Training Samples and One Test Sample');
xlabel('MDS Dimension 1'); ylabel('MDS Dimension 2');
legend('Class 1', 'Class 2', 'Class 3', 'Test Point', 'Location', 'bestoutside');
grid on; axis equal;