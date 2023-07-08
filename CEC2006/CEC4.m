function f = CEC4(P)
f(:, 2) = + 85.334407 + 0.0056858 * P(:, 2).* P(:, 5) + 0.0006262 * P(:, 1).* P(:, 4) - 0.0022053 * P(:, 3).* P(:, 5) - 92;
f(:, 3) = -85.334407 - 0.0056858 * P(:, 2).* P(:, 5) - 0.0006262 * P(:, 1).* P(:, 4) + 0.0022053 * P(:, 3).* P(:, 5);
f(:, 4) = + 80.51249 + 0.0071317 * P(:, 2).* P(:, 5) + 0.0029955 * P(:, 1).* P(:, 2) + 0.0021813 * P(:, 3).^2 - 110;
f(:, 5) = -80.51249 - 0.0071317 * P(:, 2).* P(:, 5) - 0.0029955 * P(:, 1).* P(:, 2) - 0.0021813 * P(:, 3).^2 + 90;
f(:, 6) = + 9.300961 + 0.0047026 * P(:, 3).* P(:, 5) + 0.0012547 * P(:, 1).* P(:, 3) + 0.0019085 * P(:, 3) .* P(:, 4) - 25;
f(:, 7) = -9.300961 - 0.0047026 * P(:, 3).* P(:, 5) - 0.0012547 * P(:, 1).* P(:, 3) - 0.0019085 * P(:, 3) .* P(:, 4) + 20;
f(:, 1) = 5.3578547 * P(:, 3).^2 + 0.8356891 * P(:, 1).* P(:, 5) + 37.293239 * P(:, 1) - 40792.141;