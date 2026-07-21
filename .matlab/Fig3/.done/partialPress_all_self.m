clc; clear;

% =========================================================
% F_frac sweep
% =========================================================
F_frac_all = 1e-25 : 0.003125 : 0.999999999999;
N = length(F_frac_all);

% =========================================================
% Preallocate storage
% =========================================================
a_vec = NaN(N,1);
b_vec = NaN(N,1);
c_vec = NaN(N,1);
d_vec = NaN(N,1);
e_vec = NaN(N,1);
f_vec = NaN(N,1);
g_vec = NaN(N,1);
h_vec = NaN(N,1);
ivar_vec = NaN(N,1);
jvar_vec = NaN(N,1);
P_IG_vec = NaN(N,1);

% =========================================================
% Initial guess
% =========================================================
b0 = 2e-10;       % chosen for f ~ 8e-20 at small F_frac
c0 = 1e-15;      % chosen for g ~ 1e-5 at small F_frac
x0 = [log(b0); log(c0)];

% =========================================================
% Solver options
% =========================================================
options = optimoptions('lsqnonlin', ...
    'Display','off', ...
    'MaxIterations',5000, ...
    'MaxFunctionEvaluations',10000, ...
    'StepTolerance',1e-12, ...
    'FunctionTolerance',1e-12);

% =========================================================
% Main F_frac loop
% =========================================================
for i = 1:N
    F_frac = F_frac_all(i);

    if mod(i-1,10)==0
        fprintf('Step %d/%d   F_frac = %.10g\n', i, N, F_frac);
    end

    try
        % Solve nonlinear system
        sol = lsqnonlin(@(x) eqns_log(x,F_frac), x0, [-Inf;-Inf], [Inf;Inf], options);

        % Extract solution
        B = sol(1);
        C = sol(2);
        b = exp(B);
        c = exp(C);

        % Ensure non-negative
        b = max(b,0);
        c = max(c,0);

        % Save solution for continuation
        x0 = sol;

        % -----------------------------
        % Compute derived variables (updated)
        % -----------------------------
        c_safe = max(c,1e-30);

        % Precompute core terms
        jvar = (3/2)*(b + 1.70543356735194e11*(b^2)*(c^(3/4)) + 3.56334703740337e21*(b^3)*(c^(3/2))) ...
               - 2*(c - 2.01696264109143e-20 / c_safe) ...
               + 3*(3.88741913826539e39*(b^6)*(c^3)) + 0.00425;

        ivar = 8.19244265406084e-08 * jvar / sqrt(c_safe);

        h_val = 130015309304.302 * b * c^(3/4) * sqrt(ivar);

        a_val = 9.29462654435292e13 * (b^2) * (c^(3/2));
        d_val = 1.70543356735194e11 * (b^2) * (c^(3/4));
        e_val = 3.56334703740337e21 * (b^3) * (c^(3/2));
        f_val = 3.88741913826539e39 * (b^6) * (c^3);
        g_val = 2.01696264109143e-20 / c_safe;

        % Compute P_IG
        P_IG = ((1 - F_frac)/F_frac) * ( 0.5*h_val + ivar + jvar );

        % Ensure non-negative
        a_val = max(a_val,0);
        d_val = max(d_val,0);
        e_val = max(e_val,0);
        f_val = max(f_val,0);
        g_val = max(g_val,0);
        h_val = max(h_val,0);
        ivar = max(ivar,0);
        jvar = max(jvar,0);
        P_IG = max(P_IG,0);

        % -----------------------------
        % Store results
        % -----------------------------
        a_vec(i)    = a_val;
        b_vec(i)    = b;
        c_vec(i)    = c;
        d_vec(i)    = d_val;
        e_vec(i)    = e_val;
        f_vec(i)    = f_val;
        g_vec(i)    = g_val;
        h_vec(i)    = h_val;
        ivar_vec(i) = ivar;
        jvar_vec(i) = jvar;
        P_IG_vec(i) = P_IG;

    catch
        fprintf('  ✗ Failed at F_frac = %.10g\n', F_frac);
    end
end

% =========================================================
% Store results in table
% =========================================================
Results = table(a_vec,b_vec,c_vec,d_vec,e_vec,f_vec,g_vec,h_vec, ...
                ivar_vec,jvar_vec,P_IG_vec,F_frac_all', ...
    'VariableNames',{'a','b','c','d','e','f','g','h','ivar','jvar','P_IG','F_frac'});

assignin('base','Results',Results);
fprintf('\nTable "Results" is now in the workspace.\n');

% =========================================================
% Nested function: nonlinear equations (updated)
% =========================================================
function F = eqns_log(x,F_frac)
    B = x(1);
    C = x(2);
    b = exp(B);
    c = exp(C);
    c_safe = max(c,1e-30);

    % -----------------------------
    % Precompute common terms
    % -----------------------------
    jvar = (3/2)*(b + 1.70543356735194e11*(b^2)*(c^(3/4)) + 3.56334703740337e21*(b^3)*(c^(3/2))) ...
           - 2*(c - 2.01696264109143e-20 / c_safe) ...
           + 3*(3.88741913826539e39*(b^6)*(c^3)) + 0.00425;

    j_over_sqrtc = 8.19244265406084e-08 * jvar / sqrt(c_safe);

    J = 130015309304.302 * b * c^(3/4) * sqrt(j_over_sqrtc);

    % -----------------------------
    % Equation 1
    % -----------------------------
    LHS1 = (0.000501504513540622 / F_frac) * (j_over_sqrtc + jvar);
    RHS1 = 2*( 9.29462654435292e13*b^2*c^(3/2) + 1.70543356735194e11*b^2*c^(3/4) ) ...
           + b + 3*3.56334703740337e21*b^3*c^(3/2) ...
           + 6*3.88741913826539e39*b^6*c^3 ...
           + J * (2*F_frac - 0.000501504513540622)/(2*F_frac);
    F1 = LHS1 - RHS1;

    % -----------------------------
    % Equation 2
    % -----------------------------
    F2 = 9.29462654435292e13*b^2*c^(3/2) + b + c ...
         + 1.70543356735194e11*b^2*c^(3/4) ...
         + 3.56334703740337e21*b^3*c^(3/2) ...
         + 3.88741913826539e39*b^6*c^3 ...
         + 2.01696264109143e-20 / c_safe ...
         + J + j_over_sqrtc + jvar ...
         + ((1 - F_frac)/F_frac)*(0.5*J + j_over_sqrtc + jvar) ...
         - 1.0;

    F = [F1; F2];
end