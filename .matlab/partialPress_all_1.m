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
b0 = 2e-4;       % chosen for f ~ 8e-20 at small F_frac
c0 = 1e-16;      % chosen for g ~ 1e-5 at small F_frac
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
        % Compute derived variables
        % -----------------------------
        c_safe = max(c,1e-30);

        d_val = 50567580842.629       * b^2 * c^(3/4);
        e_val = 7.01432630423488e21   * b^3 * c^(3/2);
        f_val = 1.07688142574839e40   * b^6 * c^3;
        a_val = 143355396545689       * b^2 * c^(3/2);
        g_val = 1.86489627968959e-20 / c_safe;

        jvar = (3/2)*( b + d_val + e_val ) - 2*( c - g_val ) + 3*f_val + 0.00425;
        jvar = max(jvar,0);

        ivar = 7.53238785702308e-08 * jvar / sqrt(c_safe);
        ivar = max(ivar,0);

        h_val = 162965559517.107 * b * c^(3/4) * sqrt(ivar);

        P_IG = ((1 - F_frac)/F_frac) * (0.5*h_val + ivar + jvar);

        % -----------------------------
        % Store results
        % -----------------------------
        a_vec(i)    = max(a_val,0);
        b_vec(i)    = b;
        c_vec(i)    = c;
        d_vec(i)    = max(d_val,0);
        e_vec(i)    = max(e_val,0);
        f_vec(i)    = max(f_val,0);
        g_vec(i)    = max(g_val,0);
        h_vec(i)    = max(h_val,0);
        ivar_vec(i) = ivar;
        jvar_vec(i) = jvar;
        P_IG_vec(i) = max(P_IG,0);

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
% Nested function: nonlinear equations
% =========================================================
function F = eqns_log(x,F_frac)
    B = x(1);
    C = x(2);
    b = exp(B);
    c = exp(C);
    c_safe = max(c,1e-30);

    % Power-law terms
    T1 = 50567580842.629       * b^2 * c^(3/4);
    T2 = 7.01432630423488e21   * b^3 * c^(3/2);
    T3 = 1.07688142574839e40   * b^6 * c^3;
    T4 = 143355396545689       * b^2 * c^(3/2);

    % Core j-expression
    jraw = (3/2)*( b + T1 + T2 ) - 2*( c - 1.86489627968959e-20 / c_safe ) + 3*T3 + 0.00425;
    jraw = max(jraw,1e-30);

    j_over_sqrtc = 7.53238785702308e-08 * jraw / sqrt(c_safe);
    J = 162965559517.107 * b * c_safe^(3/4) * sqrt(j_over_sqrtc);

    % -----------------------------
    % Equation 1
    % -----------------------------
    LHS1 = (0.000501504513540622 / F_frac) * (j_over_sqrtc + jraw);
    RHS1 = 2*(T4 + T1) + b + 3*T2 + 6*T3 + J * (2*F_frac - 0.000501504513540622)/(2*F_frac);
    F1 = LHS1 - RHS1;

    % -----------------------------
    % Equation 2
    % -----------------------------
    F2 = T4 + b + c + T1 + T2 + T3 + 1.86489627968959e-20 / c_safe ...
         + J + j_over_sqrtc + jraw + ((1 - F_frac)/F_frac)*(0.5*J + j_over_sqrtc + jraw) - 1.0;

    F = [F1; F2];
end