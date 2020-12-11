function A_n = pow(A,n)
    A_n = A;
    for i = 2:n
        A_n = A_n*A;
    end
end