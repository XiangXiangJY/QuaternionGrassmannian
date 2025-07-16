function complex_matrix = quaternion_to_complex(q_matrix)
    % Convert quaternion matrix to a complexified 2n x 2n matrix
    % Input:
    %   q_matrix: n x n quaternion matrix
    % Output:
    %   complex_matrix: 2n x 2n complex matrix

    [n, ~] = size(q_matrix);
    complex_matrix = zeros(2*n, 2*n); % Initialize complex matrix


    
    for i = 1:n
        for j = 1:n
            % Extract real and imaginary parts of the quaternion element
            q = q_matrix(i, j);
            r=real(q);
            a=r.w; b=r.x; c=r.y;d=r.z;% Decompose quaternion into real (a) and imaginary (b, c, d) parts


            % Fill the complexified matrix
            complex_matrix(2*i-1:2*i, 2*j-1:2*j) = ...
                [a + 1i*b, c + 1i*d; -c + 1i*d, a - 1i*b];
        end
    end
end