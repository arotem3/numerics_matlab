function [K] = Kmu(A,n)
% KMU : k means algorithm
% A : data matrix where each col is a data point
% n : number of clusters
% ---
% K : matrix of all K means

[dim,numPts] = size(A);
if dim == 2 && n < 6
    plotCheck = true;
else
    plotCheck = false;
end

if n==1
    K = mean(A,2);
    return;
end

% initialize k to be n random data points from A
temp = randperm(numPts,n); % no repeats
K = A(:,temp);

notDone = true;
while (notDone)
    K1 = K; % temp variable for checking convergence
    
    S = zeros(dim,n);
    C = zeros(1,n);
    for i=1:numPts
        a = closestK(K,n,A(:,i));
        S(:,a) = S(:,a) + A(:,i);
        C(a) = C(a) + 1;
    end
    K = S./C;
    
    if (norm(K1 - K) < 0.001) % convergence check
        notDone = false;
    end
end

if plotCheck
    colors = {'r','b','m','k','g'};
    figure
    hold on
    for i=1:numPts
        a = closestK(K,n,A(:,i));
        plot(A(1,i),A(2,i),strcat('o',colors{a}))
    end
end

end

function k = closestK(K,n,v)
min = norm(K(:,1)-v);
k = 1;
for i=2:n
    temp = norm(K(:,i)-v);
    if temp < min
        min = temp;
        k = i;
    end
end
end
