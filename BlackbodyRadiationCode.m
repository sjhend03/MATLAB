
Vp = [100, 110, 120, 130, 139.99, 149.98, 159.98, 169.98, 179.99, 189.99] ; % V
Vp_uncer = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01] ;
V = [100.047, 110.054, 120.067, 130.080, 140.082, 150.089, 160.1, 170.108, 180.123, 190.138] ; % V maybe change negative
V_uncer = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001] ;
I = [54.557, 57.754, 60.772, 63.711, 66.540, 69.318, 71.98, 74.563, 77.150, 79.610] / 1000 ; % A maybe change sign
I_uncer = [0.001, 0.001, 0.001, 0.001, 0.002, 0.001, 0.001, 0.001, 0.002, 0.001] ;
lambda_max = [738.091, 738.34, 740.585, 740.585, 739.981, 739.944, 739.302, 673.503, 571.928, 568.829] / 10^9; % m
lambda_max_uncer = [10, 10, 5, 7, 5, 6, 20, 5, 10, 5] / 10^9 ;
l = [622.717, 623.392, 621.704, 621.704, 621.704] ; % other max
I_450062 = [1010, 1550, 2285, 3100, 4400, 6000, 8000, 5200, 6700, 8400] ; 
I_450062_uncer = [70, 100, 98, 110, 100, 100, 100, 200, 200, 200] ; % counts
I_50016 = [3525, 5225, 7810, 11500, 15500, 20900, 27400, 17200, 21800, 27000] ;
I_50016_uncer = [100, 225, 200, 250, 300, 200, 300, 400, 200] ;
I_550165 = [8600, 13000, 18700, 26200, 34900, 45300, 53800, 37900, 46200, 52800] ;
I_550165_uncer = [150, 175, 300, 350, 300, 400, 300, 400, 500, 300] ;
I_60034 = [11800, 17000, 23800, 32200, 41900, 51400, 57400, 42700, 53000, 55400] ;
I_60034_uner = [300, 200, 325, 390, 300, 300, 200, 400, 300, 300] ;
it = [200, 200, 200, 200, 200, 200, 200, 100, 100, 100] ; % ms
% Resolution of 100 nm 
T_0 = 298 ;
R_0 = 167.93 ;
R_0_uncer = 0.05 ;
T = ((T_0 - 60) / R_0) .* (V ./ I) + 60 ;
T_uncer = 0 ;

P = I .* V ;
logP = log(P) ;
logT = log(T) ;
% Propagation of uncertainity
P_uncer = sqrt(I.^2.*I_uncer.^2 + V.^2.*V_uncer.^2) ;
logP_uncer = sqrt((1 ./ P.^2).*P_uncer.^2) ;
% Stefan-Boltzmann equation
% Here we are fitting ln(P) (power) vs. ln(T) (temperature)
% Our result should have a slope of 4  and an intercept of ln(epsilon*A*omega)
boltzmann_fit = wlsfit(logT, logP, logP_uncer) ;

% Wein's Law
% Here we are fitting lambda_max (highest intenisty wavelength) vs. 1 / T (Temperature)
% Our result should have a slope of 2.898*10^-3
weins_fit = wlsfit((1 ./ T), lambda_max, lambda_max_uncer) ;

% Planck's Law
% Here we are fitting ln(I) (Intensity) vs. 1 / T (Temperature)
% Our results should have a slope of -hc/lambda*k and an intercept of ln(2hc^2/epsilon)
% Function for finding chi-square and reduced chi-square with poly1 as the
% fit
function chi = findChi(x, y, y_uncer, A, B)
    chi_square = sum((y - A - B*x).^2 ./ y_uncer.^2) ;
    chi_square_reduced = chi_square / (length(y) - 2) ;
    chi = [chi_square, chi_square_reduced] ;
end
% Function for calculating linear fit by weighted least squares
function fit = wlsfit(x, y, y_uncer)
    weights = 1 ./ y_uncer.^2 ;
    delta = sum(weights)*sum(weights.*(x.^2)) - (sum(weights.*x))^2 ;
    A = (sum(weights.*x.^2)*sum(weights.*y) - sum(weights.*x)*sum(weights.*x.*y)) / delta ;
    B = (sum(weights)*sum(weights.*x.*y) - sum(weights.*x)*sum(weights.*y)) / delta ; 
    A_err = sqrt(sum(weights.*x.^2) / delta) ;
    B_err = sqrt(sum(weights) / delta) ;
    % Return in form fit = A + Bx
    fit = [A, B, A_err, B_err] ;
end

