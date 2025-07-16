% Total number of categories
num_categories = 8;

% Number of elements per category
num_elements_per_category = 10;

% Training set size (k elements)
k = 5; % Specify the number of training elements

% Number of repetitions for random sampling
m = 10; % Number of trials

% Distance matrix (assume this is already computed)
% Example: distance_matrix(i, j) represents the distance between element i and element j
% Replace with your actual distance matrix
% distance_matrix = [ ... ];

% Store accuracy for each trial
accuracy_results = zeros(m, 1);

for trial = 1:m
    % Initialize training and testing indices
    training_indices = cell(num_categories, 1);
    testing_indices = cell(num_categories, 1);

    % Randomly split the data into training and testing sets for each category
    for category = 1:num_categories
        % Get the indices for the current category
        category_indices = (category - 1) * num_elements_per_category + (1:num_elements_per_category);

        % Shuffle the indices
        shuffled_indices = category_indices(randperm(num_elements_per_category));

        % Split into training and testing sets
        training_indices{category} = shuffled_indices(1:k); % First k elements for training
        testing_indices{category} = shuffled_indices(k+1:end); % Remaining 10-k elements for testing
    end

    % Perform classification for each test element
    num_test_elements = (10 - k) * num_categories;
    predicted_labels = zeros(num_test_elements, 1); % Store predicted labels
    true_labels = repelem(1:num_categories, 10 - k)'; % True labels for test elements
    test_counter = 0; % Counter for test elements

    for category = 1:num_categories
        for test_idx = testing_indices{category}
            test_counter = test_counter + 1;

            % Shortest distances to each category's training data
            min_distances = zeros(num_categories, 1);
            for train_category = 1:num_categories
                train_indices = training_indices{train_category};
                min_distances(train_category) = min(distance_matrix(test_idx, train_indices));
            end

            % Assign the test element to the category with the minimum shortest distance
            [~, predicted_labels(test_counter)] = min(min_distances);
        end
    end

    % Compute accuracy for this trial
    accuracy_results(trial) = mean(predicted_labels == true_labels) * 100;
end

% Calculate average accuracy
average_accuracy = mean(accuracy_results);

% Display results
fprintf('Average Classification Accuracy over %d trials (Shortest Distance): %.2f%%\n', m, average_accuracy);

% Plot accuracy distribution

histogram(accuracy_results, 'Normalization', 'probability');
xlabel('Accuracy (%)');
ylabel('Frequency');
title('Accuracy Distribution over Multiple Trials (Shortest Distance)');
grid on;

% Display accuracy results
%disp('Accuracy for each trial:');
%disp(accuracy_results);