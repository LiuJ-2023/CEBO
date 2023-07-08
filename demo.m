addpath('CEC2006')

for problem_index = [1 2 4 6 7 8 9 10 18 19 24]
    for i = 1:10
        F_CEC2006{problem_index}(i,:) = CEBO(problem_index);
    end
    save Results_CEC2006
end