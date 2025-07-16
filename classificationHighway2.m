% Parameters
num_categories = 3;  % heavy, light, medium
category_sizes = [44, 165, 45];
total_samples = sum(category_sizes);
num_train_total = 192;  % Total training samples
m = 10;  % Number of trials

% Load distance matrix
a = load('highwaydistance_matrix.mat');
distance_matrix = a.distance_matrix;

% Generate ground truth labels: [1 for heavy, 2 for light, 3 for medium]
labels = [ ...
    repmat(1, category_sizes(1), 1); ...
    repmat(2, category_sizes(2), 1); ...
    repmat(3, category_sizes(3), 1)];

% Initialize accuracy results
accuracy_results = zeros(m, 1);

for trial = 1:m
    % Random split
    rand_indices = randperm(total_samples);
    train_indices = rand_indices(1:num_train_total);
    test_indices = rand_indices(num_train_total+1:end);

    train_labels = labels(train_indices);
    test_labels = labels(test_indices);
    predicted_labels = zeros(length(test_indices), 1);

    % Classify each test sample
    for t = 1:length(test_indices)
        test_idx = test_indices(t);
        avg_distances = zeros(num_categories, 1);
        for c = 1:num_categories
            class_train_idx = train_indices(train_labels == c);
            avg_distances(c) = mean(distance_matrix(test_idx, class_train_idx));
        end
        [~, predicted_labels(t)] = min(avg_distances);
    end

    % Compute trial accuracy
    accuracy_results(trial) = mean(predicted_labels == test_labels) * 100;
    fprintf('Trial %d: Accuracy = %.2f%%\n', trial, accuracy_results(trial));
end

% Final statistics
mean_accuracy = mean(accuracy_results);
std_accuracy = std(accuracy_results);

fprintf('\nAverage Accuracy over %d trials: %.2f%%\n', m, mean_accuracy);
fprintf('Standard Deviation: ±%.2f%%\n', std_accuracy);

% Plot accuracy distribution
histogram(accuracy_results, 'Normalization', 'probability');
xlabel('Accuracy (%)');
ylabel('Frequency');
title('Accuracy Distribution over Trials');
grid on;