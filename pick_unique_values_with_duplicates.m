function result = pick_half_length(D, tol)
    % Extract unique values with exactly half the length of the input vector
    % Input:
    %   D: Vector with 2n values (complex values appear with their conjugates)
    %   tol: Tolerance for numerical precision (default: 1e-6)
    % Output:
    %   result: Vector with n unique values
    
    if nargin < 2
        tol = 1e-6; % Default tolerance
    end

    % Sort values to group conjugate pairs and duplicates
    D = sort(D, 'ComparisonMethod', 'real');
    result = []; % Initialize result vector
    original_length = length(D);
    target_length = original_length / 2; % Desired length of the result vector
    
    % Step 1: Process complex numbers (remove conjugates)
    while ~isempty(D)
        current = D(1); % Take the first element
        if imag(current) ~= 0 % If it is complex
            % Find its conjugate
            conjugate_idx = find(abs(D - conj(current)) < tol, 1);
            if ~isempty(conjugate_idx)
                % Add the current value (a+bi) to result
                result = [result; current];
                % Remove both the current and its conjugate from D
                D([1, conjugate_idx]) = [];
            else
                % If no conjugate is found, treat as unique
                result = [result; current];
                D(1) = [];
            end
        else
            % If it is real, count occurrences
            real_idx = abs(D - current) < tol; % Identify all occurrences of the current real number
            count = sum(real_idx);
            
            % Keep ceil(count / 2) real values
            keep_count = ceil(count / 2);
            result = [result; repmat(current, keep_count, 1)];
            
            % Remove all occurrences of the current real number from D
            D(real_idx) = [];
        end

        % Stop processing if result reaches the desired length
        if length(result) >= target_length
            break;
        end
    end

    % Ensure the result length matches the target length
    result = result(1:target_length);
end