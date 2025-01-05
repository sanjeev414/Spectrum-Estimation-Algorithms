function opmatrix = khatrirao(A,B)    


    [m, n] = size(A);
    [p, q] = size(B);
    
    if n ~= q
        error('Matrices must have the same number of columns for Khatri-Rao product');
    end
    
    % Initialize result matrix
    opmatrix = zeros(m * p, n);
    
    % Compute the Khatri-Rao product
    for col1 = 1:n
        z = [];
        for i = 1:m
            for j = 1:p
                z = [z;A(i,col1)*B(j,col1)];

            end
        end
        opmatrix(:,col1) = z;
    end

end