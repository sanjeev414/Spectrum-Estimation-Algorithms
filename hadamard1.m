function result = hadamard1(A,B)
    % Multiply two matrices where each element of A is multiplied with all
    [M,N] = size(A);
    [P,Q] = size(B);
    
    % Initialize the resulting matrix
    result = zeros(M * P, N * Q);

    % Compute the custom product
    for i = 1:M
        for j = 1:N
            % Compute the submatrix for element A(i, j)
            submatrix = A(i, j) * B;

            % Determine placement in the result matrix
            row_start = (i - 1) * P + 1;
            row_end = row_start + P - 1;
            col_start = (j - 1) * Q + 1;
            col_end = col_start + Q - 1;

            % Place the submatrix
            result((row_start:row_end), (col_start:col_end)) = submatrix;
        end
    end

end
