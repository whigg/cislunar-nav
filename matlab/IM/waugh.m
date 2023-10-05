         0.779509185574721
       -0.0148441669777531
         0.845719600164455


             % perform a cubic-interpolation backtracking line search
    phi0 = F(x);
    dphi0 = dot(dF, p);
    phi_prev = 0;
    alpha_prev = 0;
    for j=1:max_iter
        phi = F(x + alpha * p);

        if phi <= phi0 + c1 * alpha * dphi0
            break
        else
            if j == 1
                b = dphi0;
                a = (phi - dphi0 * alpha - phi0) / alpha^2;
                alpha_new = -b / (2*a);
            else
                A = [alpha^3 alpha^2; alpha_prev^3 alpha_prev^2];
                B = [phi - dphi0 * alpha - phi0; phi_prev - dphi0 * alpha_prev - phi0];
                sol = A \ B;
                a = sol(1); b = sol(2);
                alpha_new = (-b + sqrt(b^2 - 3*a*dphi0))/(3*a);
            end

            alpha_prev = alpha;
            phi_prev = phi;

            if alpha_new > 0.95 * alpha
                alpha = 0.95 * alpha;
            elseif alpha_new < 0.05 * alpha
                alpha = 0.05 * alpha;
            else
                alpha = alpha_new;
            end
        end

        if j == max_iter
            throw(MException('Optimizer:BacktrackingLineSearch', ...
                'Failed to converge to a in %d iterations.', j));
        end
    end