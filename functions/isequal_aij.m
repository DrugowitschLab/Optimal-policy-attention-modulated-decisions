function output = isequal_aij(mat1,mat2,flag)
%isequal, but takes care of precision error
% flag
% 1: give binary output
% 2: give output for every element
approxZero = 1.0e-8;
temp=abs(mat1-mat2);
if flag==1
    if any(temp(:) > approxZero)
        output = false;
    else
        output = true;
    end
elseif flag==2
    output = temp < approxZero;
end

end

