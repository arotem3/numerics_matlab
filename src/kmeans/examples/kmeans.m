n = 20;
A = [rand(n,2) + 1; rand(n,2) - 1]; %[n,2]
plot(A(:,1),A(:,2),'o')

notDone = true;
k1 = A(randi([1,2*n]),:); % initialize k1,k2
k2 = A(randi([1,2*n]),:);

t = linspace(-3,3,10);
v = f(k1,k2,t);
hold on
plot(k1(1),k1(2),'rs')
plot(k2(1),k2(2),'rs')
plot(v(:,1),v(:,2))
xlim([-1,2])
ylim([-1,2])
pause(1)
while(notDone)
    k11 = k1;
    k21 = k2;
    
    temp1 = [0,0];
    count1 = 0;
    temp2 = [0,0];
    count2 = 0;
    for i=1:2*n
        if leq(k1,k2,A(i,:))
            temp1 = temp1 + A(i,:);
            count1 = count1 + 1;
        else
            temp2 = temp2 + A(i,:);
            count2 = count2 + 1;
        end
    end
    k1 = temp1/count1;
    k2 = temp2/count2;
    if (norm(k1-k11) < 0.01 && norm(k2-k21) < 0.01)
        notDone = false;
    end
    
    clf;
    plot(A(:,1),A(:,2),'o')
    hold on
    v = f(k1,k2,t);
    plot(k1(1),k1(2),'rs')
    plot(k2(1),k2(2),'rs')
    plot(v(:,1),v(:,2))
    xlim([-1,2])
    ylim([-1,2])
    pause(1)
end

function v = f(x1,x2,t)
v(:,1) = -(x1(2) - x2(2)) * t + (x1(1) + x2(1))/2;
v(:,2) = (x1(1) - x2(1)) * t + (x1(2) + x2(2))/2;
end

function z = leq(k1,k2,x)
y = (k1(2) - k2(2))/(k2(1) - k1(1)) * (x(1) - (k1(1) + k2(1))/2) + (k1(2) + k2(2))/2;
if x(2) <= y, z = true;
else, z = false;
end
end