function f = Fitness_CEC2006(x,index)
if index == 1
    f = CEC1(x);
elseif index == 2
    f = CEC2(x);
elseif index == 4
    f = CEC4(x);
elseif index == 6
    f = CEC6(x);
elseif index == 7
    f = CEC7(x);
elseif index == 8
    f = CEC8(x);
elseif index == 9
    f = CEC9(x);
elseif index == 10
    f = CEC10(x);
elseif index == 12
    f = CEC12(x);
elseif index == 18
    f = CEC18(x);
elseif index == 19
    f = CEC19(x);
elseif index == 24
    f = CEC24(x);
end