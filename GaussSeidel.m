% GaussSeidel function
% Kevin(Yinuo) Huang CID:01051134 16:13-18:20 29/03/2016
% Ax=B

function x=GaussSeidel(A,B,tolerence,MaxNumOfIter)
    % A and B are input matrix, err is the tolerence, NumOfIter is the
    % number of iterations
    % x is the output solution
    % take the size of two matrix
    [n,m]=size(A);
    [u,~]=size(B);
    

    counter=1; % actual number of iterations used
    x0=zeros(n,1); % we need a zero matrix to build in
    x=zeros(n,1); % initialise output matrix x

    
    while (counter<MaxNumOfIter) % loop ends when exceed max no. of iterations 
        %Gauss-Seidel Iteration
        for i=1:n
            I=[1:i-1 i+1:n];
            x(i)=(B(i)-A(i,I)*x(I))/A(i,i);
        end
        %calculate error and compare with tolerence entered
        esp=(abs(x(i)-x0)/abs(x(i)));
        %break if good enough
        if max(esp)<tolerence
           break;
        end
        x0=x; % assign the new x to x0
        counter=counter+1;% no. of iterations
    end

    
