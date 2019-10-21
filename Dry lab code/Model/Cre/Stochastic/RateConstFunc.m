function RateConst = RateConstFunc(x, t, c)
% Edit the function according to the need
    v = length(c);                                                  % Number of reactions
    RateConst = zeros(1, v);
    
% Edit Rate parameters according to the need
    RateConst(1) = c(1)*x(1)*x(14);
    RateConst(2) = c(2)*x(3);
    RateConst(3) = c(3)*x(2)*x(14);
    RateConst(4) = c(4)*x(4);
    RateConst(5) = c(5)*x(3)*x(14);
    RateConst(6) = c(6)*x(5);
    RateConst(7) = c(7)*x(4)*x(14);
    RateConst(8) = c(8)*x(6);
    RateConst(9) = c(9)*x(5)*x(6);
    RateConst(10) = c(10)*x(7);
    RateConst(11) = c(11)*x(7);
    RateConst(12) = c(12)*x(8)*x(9);
    RateConst(13) = c(13)*x(8);
    RateConst(14) = c(14)*x(10)*x(14);
    RateConst(15) = c(15)*x(9);
    RateConst(16) = c(16)*x(11)*x(14);
    RateConst(17) = c(17)*x(10);
    RateConst(18) = c(18)*x(12)*x(14);
    RateConst(19) = c(19)*x(11);
    RateConst(20) = c(20)*x(13)*x(14);
    RateConst(21) = c(21)*x(1);
    RateConst(22) = c(22)*x(15);
end