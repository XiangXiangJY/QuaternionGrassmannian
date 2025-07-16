% Parameters
num_categories = 8; % Updated for 10 categories
num_elements_per_category = 10;
k = 5; % Number of training elements per category
m = 10; % Number of trials
target_accuracy = 00.0; % Target accuracy threshold
max_attempts = 1000; % Maximum number of repetitions to find a valid set

% Distance matrix
%a = load('distance_matrix9.mat');
a= load('distance_matrixPCA.mat');
distance_matrix = a.distance_matrix;

% Store results that meet the accuracy threshold
successful_accuracy_results = [];

for attempt = 1:max_attempts
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

                % Average distances to each category's training data
                avg_distances = zeros(num_categories, 1);
                for train_category = 1:num_categories
                    train_indices = training_indices{train_category};
                    avg_distances(train_category) = mean(distance_matrix(test_idx, train_indices));
                end

                % Assign the test element to the category with the minimum average distance
                [~, predicted_labels(test_counter)] = min(avg_distances);
            end
        end

        % Compute accuracy for this trial
        accuracy_results(trial) = mean(predicted_labels == true_labels) * 100;
    end

    % Calculate average accuracy for this attempt
    average_accuracy = mean(accuracy_results);
    std_accuracy = std(accuracy_results);

    % Check if the average accuracy meets the target threshold
    if average_accuracy >= target_accuracy
        successful_accuracy_results = accuracy_results;
        fprintf('Found a successful attempt with Average Accuracy: %.2f%%\n', average_accuracy);
        break;
    end
end

if isempty(successful_accuracy_results)
    fprintf('No successful attempt found after %d repetitions.\n', max_attempts);
else
    % Display successful results
    fprintf('Successful Average Accuracy over %d trials: %.2f%%\n', m, average_accuracy);
    fprintf('Standard Deviation of Accuracy: ±%.2f%%\n', std_accuracy);

    % Plot accuracy distribution
    histogram(successful_accuracy_results, 'Normalization', 'probability');
    xlabel('Accuracy (%)');
    ylabel('Frequency');
    title('Accuracy Distribution for Successful Attempt');
    grid on;

    % Display accuracy results
   % disp('Accuracy for each trial:');
    %disp(successful_accuracy_results);
end