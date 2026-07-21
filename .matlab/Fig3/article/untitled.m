if GR == 0.01
    b = exp(xsol(1));
    c = exp(xsol(2));
    g = 1.86489627968959e-20 / max(c,1e-30);

    F = eqns_bc_new(xsol,F_frac,GR);

    fprintf('\n======= DIAGNOSTIC @ GR=0.01 =======\n');
    fprintf('b = %.16e\n', b);
    fprintf('c = %.16e\n', c);
    fprintf('g = %.16e\n', g);
    fprintf('Residual eq1 = %.3e\n', F(1));
    fprintf('Residual eq2 = %.3e\n', F(2));
    fprintf('exitflag = %d\n', exitflag(i));
    fprintf('====================================\n\n');
end