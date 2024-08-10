horner("input3.txt")
newtonRunner("input3.txt")
cramer("input.txt")


%the next three functions behave similarly, they are meant to streamline
%testing. cramer(), newtonRunner(), and horner() all take the path and file
%name of the desired input as parameters, and will display the
%corresponding results
function cramer(fileName)
    [A, b] = cramerInput(fileName);
    cramerDriver(A, b);
end

function newtonRunner(fileName)
    [x, y, e, N] = newtonInput(fileName);
    newtonDriver(x, y, e, N);
end

function horner(fileName)
    [x,y] = hornerInput(fileName);
    p = hornerDriver(x,y);
    disp("P(" + num2str(y)+") = " + num2str(p(1)))
    for i = 2:size(p)
        disp("P^(" + num2str(i-1) + ") (" + num2str(y) + ") = " +num2str(p(i)))
    end
end
%%%


%handles calling gaussian() and determinant() for the different matrices,
%and uses those values to calculate the solution via Cramer's rule
%handles displaying everything
%A and b - the matrix and vector
function cramerDriver(A, b)
    copyOfA = A;
    [m, rs] = gaussian(A);
    det = determinant(m, rs);
    dets = zeros(size(A,1));
    for i = 1:size(A,1)
        A = copyOfA;
        A(:,i) = b;
        [m, rs] = gaussian(A);
        dets(i) = determinant(m,rs);
    end
    error = false;
    for i=1:size(dets)
        if(dets(i) == false)
            error = true;
        end
    end
    if(error == true)
        disp("Error");
    else
        disp("Determinant A = " + num2str(det))
        for i=1:size(dets)
            disp("Determinant A" + num2str(i) +" = " + num2str(dets(i)));
        end
        for i = 1:size(dets)
            disp("x" + num2str(i) + " = " + num2str(dets(i)/det));
        end
    end
end

          




%takes an upper triangular matrix and the amount of row swaps needed and
%computes the determinant
%A - upper triangular matrix
%rowSwaps - the amount of times row swaps were performed
%returns the determinant, or false if impossible
function det = determinant(A, rowSwaps)
    if(rowSwaps == -1)
        det = false;
    else
        det = 1;
        for i=1:size(A,1)
            det = det * A(i,i);
            det = det * (-1)^rowSwaps;
        end
    end
end

%performs gaussian elimination on A to get it upper triangular
%A - matrix to row reduce
%returns the upper triangular matrix and the amount of row swaps needed
function [A, rowSwaps] = gaussian(A)
    rowSwaps = 0;
    for i = 1:size(A,1)
        if(A(i,i) == 0)
            rowToSwapIdx = idxToRowSwap(A, i);
            if(rowToSwapIdx < 0)
                rowSwaps = -1;
                break;
            end
            A([i rowToSwapIdx],:) = A([rowToSwapIdx i],:);
            rowSwaps = rowSwaps + 1;
        end
        for j = (i + 1):size(A,1)
            scalar = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - scalar * A(i,:);
        end
        
    end
end

%used for swapping rows in Gaussian Elimination
%A - the matrix in question
%row - the index of the row needing swapping
%returns the index of the row to swap with, or -1 if no such row exists
function idx = idxToRowSwap(A, row)
    idx = -1;
    for i = row+1:size(A,1)
        if(A(i, row) ~= 0)
            idx = i;
            break;
        end
    end
end

%reads input files for Cramer's method
%fileName - str with the path and file name of the input
%returns A and b
function [A, b] = cramerInput(fileName)
    fid = fopen(fileName, "r");
    input = fscanf(fid, '%f');
    matSize = input(1);
    input(1) = [];
    A = zeros(matSize, matSize);
    b = zeros(matSize, 1);
    for i = 1:matSize
        for j = 1:matSize
            A(i,j) = input(matSize*(i-1)+ j);
        end
    end
    for i = 1:matSize
        b(i) = input(end);
        input(end) = [];
    end
    b = flip(b);
end

%runs and processes newton, displaying the results
%A - coefficient vector
%x0 - initial guess
%e - error bound
%N - maximum iterations
function newtonDriver(A, x0, e, N)
    r = newton(A, x0, e, N);
    if(r == false)
        disp("did not converge")
    else
        disp("root = " + num2str(r))
    end
end

%utilizes horner's to approximate the root
%A - coefficient vector
%x0 - initial guess
%e - error bound
%N - maximum iterations
%returns x1 AKA the approximation, or false if N is exceeded
function r = newton(A, x0, e, N)
    x1 = 0;
    for i = 1:N
        hornerAnswers = hornerDriver(A, x0);
        x1 = x0 - (hornerAnswers(1)/hornerAnswers(2));
        if(abs(x1 - x0) <= e)
            break;
        end
        x0 = x1;
    end
    if(abs(x1 - x0) <= e)
        r = x1;
    else
        r = false;
    end
end

%fileName - the path and file name of the input
%x - the coefficients, from the greatest degree, descending
%y - x0
function [x,y, e, N] = newtonInput(fileName)
    fid = fopen(fileName, 'r');
    A = fscanf(fid, '%f');
    N = A(end);
    A(end) = [];
    e = A(end);
    A(end) = [];
    y = A(end);
    A(end) = [];
    A(1) = [];
    x = flip(A);
end

%uses synthetic division to run Horner's algo
%A - vector representing coefficients
%x0 - x value to solve at
%returns a vector of the multiple derivatives at x0, with the first element
%being the solution
function x = hornerDriver(A, x0)
    coeffs = A;
    
    x = zeros(size(A));
    for i = 1:size(A)
        [synRet, p] = syntheticDivision(coeffs, x0);
        if(i > 1)
            p = p * factorial(i-1);
        end
        x(i) = p;
        
        coeffs = synRet;
    end
end




%A - vector containing polynomials, starting with the highest order
%x0 - initial guess
%x - answer vector
%y - return value
function [x,y] = syntheticDivision(A, x0)
    ret = zeros(size(A));
    ret(1)= A(1);
    for i = 2:size(A)
        toAdd = ret(i-1) * x0;
        ret(i) = A(i) + toAdd;
    end
    y = ret(end);
    ret(end) = [];
    x = ret;
end

%fileName - the path and file name of the input
%x - the coefficients, from the greatest degree, descending
%y - x0
function [x,y] = hornerInput(fileName)
    fid = fopen(fileName, 'r');
    A = fscanf(fid, '%f');
    y = A(end);
    A(end) = [];
    A(1) = [];
    x = flip(A);
end
