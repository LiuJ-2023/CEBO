function f = CEC9(P)
f(:, 2) = -127 + 2 * P(:, 1).^2 + 3 * P(:, 2).^4 + P(:, 3) + 4 * P(:, 4).^2 + 5 * P(:, 5);
f(:, 3) = -282 + 7 * P(:, 1) + 3 * P(:, 2) + 10 * P(:, 3).^2 + P(:, 4) - P(:, 5);
f(:, 4) = -196 + 23 * P(:, 1) + P(:, 2).^2 + 6 * P(:, 6).^2 - 8 * P(:, 7);
f(:, 5) = 4 * P(:, 1).^2 + P(:, 2).^2 - 3 * P(:, 1).* P(:, 2) + 2 * P(:, 3).^2 + 5 * P(:, 6) - 11 * P(:, 7);
f(:, 1) = (P(:, 1) - 10).^2 + 5 * (P(:, 2) - 12).^2 + P(:, 3).^4 + 3 * (P(:, 4) - 11).^2 + 10 * P(:, 5).^6 + ...
    7 * P(:, 6).^2 + P(:, 7).^4 - 4 * P(:, 6).* P(:, 7) - 10 * P(:, 6) - 8 * P(:, 7);