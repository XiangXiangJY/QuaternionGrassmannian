function eigenvalues = quaternion_eigenvalues(q_matrix)
    % Compute the eigenvalues of a quaternion matrix using complexified representation
    % Input:
    %   q_matrix: n x n quaternion matrix
    % Output:
    %   eigenvalues: Eigenvalues of the quaternion matrix

    % Convert quaternion matrix to complexified representation
    complex_matrix = quaternion_to_complex(q_matrix);

    % Compute eigenvalues of the complex matrix
    complex_eigenvalues = eig(complex_matrix);

    % Extract real and imaginary parts for interpretation
    eigenvalues = complex_eigenvalues; % Eigenvalues in the complexified form
end

