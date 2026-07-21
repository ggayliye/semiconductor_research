clc;
clear;

%% =====================================================
% FIXED PARAMETERS
%% =====================================================
F_frac_init = 1e-7;       % very small fraction
GR_vec = [0.01:0.001125:0.1, ...   % stop before 0.1
          0.11125:0.01125:1, ...   % stop before 1.0
          1.1125:0.1125:10, ...    % stop before 10
          11.125:1.125:100];       % stop before 100

N = length(GR_vec);

%% =====================================================
% SOLVER OPTIONS (Levenberg–Marquardt)
%% =====================================================
options = optimoptions('lsqnonlin', ...
    'Algorithm','levenberg-marquardt', ...
    'Display','off', ...
    'MaxIterations',3000, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12);

%% =====================================================
% STORAGE
%% =====================================================
vars = {'a','b','c','d','e','f','g','h','ivar','jvar','P_IG','exitflag'};
for v = vars
    eval([v{1}, '_vec = NaN(N,1);']);
end

%% =====================================================
% INITIAL GUESS (log-space)
%% =====================================================
x0 = log([7.49903e-3; 2e-13]);   % b0 and c0

%% =====================================================
% MAIN GR LOOP
%% =====================================================
F_frac = F_frac_init;

for i = 1:N
    GR = GR_vec(i);
    fprintf('GR step %d/%d | GR = %.6g\n', i, N, GR);

    try
        [xsol,~,~,ef] = lsqnonlin(@(x) eqns_bc_safe_dynamic(x,F_frac,GR), x0, [], [], options);
        exitflag_vec(i) = ef;
    catch ME
        fprintf('  ✗ Solver failed at step %d (GR=%.6g): %s\n', i, GR, ME.message);
        continue
    end

    %% ---- Extract unknowns ----
    b = exp(xsol(1));
    c = exp(xsol(2));
    c_safe = max(c,1e-50);

    %% ---- Compute derived variables ----
    a = 143355396545689*(b^2)*(c^(3/2));
    d = 50567580842.629*(b^2)*(c^(3/4));
    e = 7.01432630423488E+21*(b^3)*(c^(3/2));
    f = 1.07688142574839E+40*(b^6)*(c^3);
    ivar = 0.0;  % ivar is now always zero
    h = 162965559517.107*b*(c^(3/4)) * sqrt(max(ivar,0) + 1e-40);
    jvar = (3/2)*(b + d + e) - 2*(c - 1.86489627968959E-20/c_safe) + 3*f ...
           - (3/2)*(0.0005 - (4/3)*GR*0.0005/2);
    g = 1.86489627968959E-20 / c_safe;
    P_IG = ((1 - F_frac)/F_frac) * (0.5*h + ivar + jvar);

    %% ---- Check for negative values ----
    if any([a,d,e,f,g,h,jvar,P_IG]<0)
        fprintf('  ⚠ Warning: Negative value detected at GR=%.6g\n', GR);
    end

    %% ---- Store ----
    a_vec(i)=a; b_vec(i)=b; c_vec(i)=c; d_vec(i)=d; e_vec(i)=e;
    f_vec(i)=f; g_vec(i)=g; h_vec(i)=h; ivar_vec(i)=ivar;
    jvar_vec(i)=jvar; P_IG_vec(i)=P_IG;

    %% ---- Continuation ----
    x0 = xsol;
end

%% =====================================================
% RESULTS TABLE
%% =====================================================
Results = table(GR_vec', a_vec, b_vec, c_vec, d_vec, e_vec, f_vec, g_vec, ...
                h_vec, ivar_vec, jvar_vec, P_IG_vec, exitflag_vec, ...
                'VariableNames',{'GR','a','b','c','d','e','f','g','h','ivar','jvar','P_IG','exitflag'});

disp('✅ CLEAN GR SWEEP COMPLETE');

%% =====================================================
% NONLINEAR EQUATIONS FUNCTION
%% =====================================================
function F = eqns_bc_safe_dynamic(x,F_frac,GR)
b = exp(x(1));
c = exp(x(2));
c_safe = max(c,1e-50);

% Derived variables
a = 143355396545689*(b^2)*(c^(3/2));
d = 50567580842.629*(b^2)*(c^(3/4));
e = 7.01432630423488E+21*(b^3)*(c^(3/2));
f = 1.07688142574839E+40*(b^6)*(c^3);
ivar = 0.0;  % ivar is now always zero
h = 162965559517.107*b*(c^(3/4)) * sqrt(max(ivar,0) + 1e-40);
jvar = (3/2)*(b + d + e) - 2*(c - 1.86489627968959E-20/c_safe) + 3*f ...
       - (3/2)*(0.0005 - (4/3)*GR*0.0005/2);
g = 1.86489627968959E-20 / c_safe;
P_IG = ((1 - F_frac)/F_frac) * (0.5*h + ivar + jvar);

% Nonlinear equations
prefactor_nom = 0.0005/(0.9995 - (GR*0.0005)/2);
prefactor = prefactor_nom / F_frac;
RHS_factor = (2*F_frac - prefactor_nom)/(2*F_frac);

eq1 = prefactor*(ivar + jvar) - (2*(a + d) + b + 3*e + 6*f + h*RHS_factor);
eq2 = a + b + c + d + e + f + g + h + ivar + jvar + ((1-F_frac)/F_frac)*(0.5*h + ivar + jvar) - 1.0;

F = [eq1; eq2];
end