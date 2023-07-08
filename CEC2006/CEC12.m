function f = CEC12(P)
[popsize,~] = size(P);
f(:,1) = -(100 - (P(:, 1) - 5).^2 - (P(:, 2) - 5).^2 - (P(:, 3) - 5).^2)/100;
for j = 1:popsize

    f(j+1, 1) = min(sum((repmat(P(j, :), 9 * 9 * 9, 1) - aaa).^2, 2)) - 0.0625;

end