function dv = orbitandlyapdyn(t,v,f,dfdx,Q)
    %ORBITANDLYAPDYN Describes the joint state dynamics and continuous-time
    %Lyapunov equations. For use with ODE45 or other propagator.
    %
    %   Inputs: (dims)
    %    - t; time step
    %    - v; vectorized full state
    %    - f; dynamics function
    %    - dfdx; dynamics partials function
    %    - Q; (nxn) process noise matrix (may also take form LQL')

    n = size(Q, 1);
    x_ = v(1:n);
    P_ = v(n+1:end);
    A  = dfdx(x_);
    dv = [f(t, x_); lyapunov(P_, A, Q)];
end
