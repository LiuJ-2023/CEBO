function f = CEC10(P)
f(:, 2) = -1 + 0.0025 * (P(:, 4) + P(:, 6));
f(:, 3) = -1 + 0.0025 * (P(:, 5) + P(:, 7) - P(:, 4));
f(:, 4) = -1 + 0.01 * (P(:, 8) - P(:, 5));
f(:, 5) = -P(:, 1).* P(:, 6) + 833.33252 * P(:, 4) + 100 * P(:, 1) - 83333.333;
f(:, 6) = -P(:, 2).* P(:, 7) + 1250 * P(:, 5) + P(:, 2).* P(:, 4) - 1250 * P(:, 4);
f(:, 7) = -P(:, 3).* P(:, 8) + 1250000 + P(:, 3).* P(:, 5) - 2500 * P(:, 5);
f(:, 1) = P(:, 1) + P(:, 2) + P(:, 3);