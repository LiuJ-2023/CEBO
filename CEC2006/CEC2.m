function f = CEC2(P)
[popsize,~] = size(P);
f(:, 2) = 0.75 - prod(P, 2);
f(:, 3) = sum(P')' - 7.5 * size(P, 2);
f(:, 1) = -abs(sum((cos(P).^4), 2) - 2 * prod((cos(P).^2), 2)) ./ sqrt(1E-30 + sum(repmat(1 : size(P, 2), popsize, 1) .* (P.^2), 2));