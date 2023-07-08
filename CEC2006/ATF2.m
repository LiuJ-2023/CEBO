function f = ATF2(P)
[m,n] = size(P);
for i = 1:m
    f(i, 2) = min(3*(P(i,1) + 7)^2 + P(i,2)^2 - 0.4,(P(i,1) + 8)^2 + (P(i,2) - 3)^2 - 2);
    f(i, 1) = P(i,1)^2 + P(i,2)^2;% - 10*cos(2*pi*P(i,1)) -  10*cos(2*pi*P(i,2)) + 10;
end