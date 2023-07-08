function f = CEC1(P)
f(:, 2) = 2 * P(:, 1) + 2 * P(:, 2) + P(:, 10) + P(:, 11) - 10;
f(:, 3) = 2 * P(:, 1) + 2 * P(:, 3) + P(:, 10) + P(:, 12) - 10;
f(:, 4) = 2 * P(:, 2) + 2 * P(:, 3) + P(:, 11) + P(:, 12) - 10;
f(:, 5) = -8 * P(:, 1) + P(:, 10);
f(:, 6) = -8 * P(:, 2) + P(:, 11);
f(:, 7) = -8 * P(:, 3) + P(:, 12);
f(:, 8) = -2 * P(:, 4) - P(:, 5) + P(:, 10);
f(:, 9) = -2 * P(:, 6) - P(:, 7) + P(:, 11);
f(:, 10) = -2 * P(:, 8) - P(:, 9) + P(:, 12);
f(:, 1) = 5 * sum(P(:, 1 : 4), 2) - 5 * sum(P(:, 1 : 4).^2, 2) - sum(P(:, 5 : 13), 2);