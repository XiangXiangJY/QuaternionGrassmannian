function Z = QuaternionPCA_Column(Q, t)
% QuaternionPCA_Column - Perform PCA on quaternion column vectors
% Q : m × N quaternion matrix (each column is a quaternion vector of length m)
% t : reduced dimension (t < m)
% Z : t × N reduced quaternion vectors
% U : m × t quaternion projection matrix
Q=quaternion([1,2;3,4],[1,2;3,4],[1,2;3,4],[1,2;3,4])
    [m, N] = size(Q);
t=1;
    % Step 1: Mean Centering
    Q_mean = mean(Q, 2);                          % m × 1
    Q_c = Q - Q_mean * ones(1, N);               % m × N centered quaternion matrix

    % Step 2: Quaternion covariance matrix: C = Q_c * Q_c^H
    Q_cR=QuaternionReal(Q_c);
    Q_cHR=transQ(Q_cR);
                                % Hermitian transpose (conjugate + transpose)
    C = timesQ(Q_cHR,Q_cR);                               % m × m quaternion covariance matrix

    % Step 3: Convert quaternion covariance to complex 2m × 2m matrix
    [U,D,V]= svdQ(C);       % complex surrogate matrix

       % sort by eigenvalue
       VQ=RealQuaternion(V);
    V_t = VQ(:, 1:t);                 % select top t components

    % Step 5: Convert eigenvectors back to quaternion form
   V_tR=QuaternionReal(V_t);
    % Step 6: Project centered quaternion vectors to reduced space
    ZR = timesQ(Q_cR,V_tR);
  Z=RealQuaternion(ZR);% t × N reduced quaternion matrix
end