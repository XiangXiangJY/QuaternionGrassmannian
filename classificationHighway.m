% Parameters
num_categories = 3;  % heavy, light, medium
category_sizes = [44, 165, 45];
total_samples = sum(category_sizes);
num_train_total = 192;  % Total number of training samples
m = 10;  % Number of trials
target_accuracy = 0.0;  % No threshold; collect results
max_attempts = 1000;

% Load distance matrix
a = load('highwaydistance_matrix.mat');
distance_matrix = a.distance_matrix;

% Generate ground truth labels in order: heavy, light, medium
labels = [ ...
    repmat(1, category_sizes(1), 1); ...
    repmat(2, category_sizes(2), 1); ...
    repmat(3, category_sizes(3), 1)];

% Store successful results
successful_accuracy_results = [];

for attempt = 1:max_attempts
    accuracy_results = zeros(m, 1);

    for trial = 1:m
        % Randomly select training indices
        rand_indices = randperm(total_samples);
        train_indices = rand_indices(1:num_train_total);
        test_indices = rand_indices(num_train_total+1:end);

        train_labels = labels(train_indices);
        test_labels = labels(test_indices);
        predicted_labels = zeros(length(test_indices), 1);

        % Classification by nearest average distance
        for t = 1:length(test_indices)
            test_idx = test_indices(t);
            avg_distances = zeros(num_categories, 1);
            for c = 1:num_categories
                class_train_idx = train_indices(train_labels == c);
                avg_distances(c) = mean(distance_matrix(test_idx, class_train_idx));
            end
            [~, predicted_labels(t)] = min(avg_distances);
        end

        % Accuracy for this trial
        accuracy_results(trial) = mean(predicted_labels == test_labels) * 100;
    end

    % Evaluate attempt
    average_accuracy = mean(accuracy_results);
    std_accuracy = std(accuracy_results);

    if average_accuracy >= target_accuracy
        successful_accuracy_results = accuracy_results;
        fprintf('Found a successful attempt with Average Accuracy: %.2f%%\n', average_accuracy);
        break;
    end
end

% Display results
if isempty(successful_accuracy_results)
    fprintf('No successful attempt found after %d repetitions.\n', max_attempts);
else
    fprintf('Successful Average Accuracy over %d trials: %.2f%%\n', m, average_accuracy);
    fprintf('Standard Deviation of Accuracy: ±%.2f%%\n', std_accuracy);

    % Plot accuracy distribution
    histogram(successful_accuracy_results, 'Normalization', 'probability');
    xlabel('Accuracy (%)');
    ylabel('Frequency');
    title('Accuracy Distribution for Successful Attempt');
    grid on;
end