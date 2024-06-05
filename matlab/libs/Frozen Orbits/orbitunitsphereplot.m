function orbitunitsphereplot(oes,t0,moon,earth,frame)
%ORBITUNITSPHEREPLOT Plots orbits on a unit sphere to qualitatively analyze
%geometry for dilution of precision.

figure();
[xx, yy, zz] = ellipsoid(0, 0, 0, 1, 1, 1, 15);
globe = surf(xx, yy, zz);
set(globe, 'FaceColor', 'none', 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.5);
hold on;

opts = odeset("RelTol", 1e-8, "AbsTol", 1e-9);
R = rotx(-6.76*pi/180);
for oe=oes
    T = 2*pi * sqrt(oe.a^3/moon.GM);
    [r,v] = oe2rv(oe.a,oe.e,oe.i,oe.RAAN,oe.w,oe.f,moon.GM);
    x0 = cspice_sxform('MOON_OP', 'J2000', t0) * [r; v];
    ts = t0:600:t0+T;
    [~,X] = ode45(@(t,x) orbitaldynamics(t,x,moon,50,earth), ts, x0, opts);
    X = X';
    for i=1:length(ts)
        X(:,i) = cspice_sxform('J2000', frame, ts(i)) * X(:,i);
    end
    X = X(1:3,:);
    X = R*X ./ sqrt(sum(X.^2,1));
    X = X(:,X(3,:) < 0);
    scatter3(X(1,:), X(2,:), X(3,:),'filled', 'o');
end
hold off; axis equal;
end

