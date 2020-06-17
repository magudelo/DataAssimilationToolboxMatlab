function sqrtm = sqrtmSchur(matrix)
            
    if isequal(matrix,diag(diag(matrix)))       %diagonal matrix
        
        sqrtm = diag(sqrt(diag(matrix)));
        
    else
        n = length(matrix);
        [U, T] = schur(matrix,'complex');  %matrix=UTU'

        if any(diag(T) <= eps)
            warning('DA:dautils:sqrtmSchur:notPosDef','Matrix is numerically singular.');
        end
            
    	%R= square root of T
        if isequal(T,diag(diag(T)))      % Check if T is diagonal.
    
            R = diag(sqrt(diag(T)));     % Square root always exists. 
    
        else
            % Compute upper triangular square root R of T, a column at a time.
            R = zeros(n);  
            for j=1:n
                R(j,j) = sqrt(T(j,j));
                for i=j-1:-1:1
                    k = i+1:j-1;
                    s = R(i,k)*R(k,j);
                    R(i,j) = (T(i,j) - s)/(R(i,i) + R(j,j));
                end
            end
        end
        
        sqrtm = U * R * U';
    end
            
end    