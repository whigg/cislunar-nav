function dy = dynamicsSTM_CR3BP(t, y, mu)
%This file contains the right hand side for the system which defines the
%state transition matrix for the simple unperturbed CR3BP dynamics. So "dy"
%is a vector with 42 components
%-------------------------------------------------------------------------
%The structure of the system is
% STM’ = A*STM, x’ = f(x)
%
%Where STM (here called: A) are 12X12 matrices and STM’ is the derivative of STM.
%A is the derivative of the N-body vector field "f" and has the form
%
%      |  0     I  |
% A =  | A21   A22 |
%
%where all submatrices are 3X3.
%The non-zero matricies are "I" = Identity Matrix or constant matricies
%which can be found in the report.
%All of this is combined onto one system of 42 1st order ODEs. The right
%hand side for the system is coded in the remainder of the file.
%
%Author: Luigi De Maria, 2022
%(readaptation of original 2BP code)
%------------------------------------------------------------------------
s = y(7:42); % STM
x = y(1:6); % Sates
%-------------------------------------------------------------------------
% Now compute the State Matrix A:
%Zero and identity Matricies Allocation
O = zeros(3);
I = eye(3);
%%%%%%%%%%%%%%%%%%%%%%% 1ST DERIVATIVES OF g(r) BY r %%%%%%%%%%%%%%%%%%%%%%
%---------------------------------- g1 ------------------------------------
dg1dx = @(x,y,z) 1 - (1-mu) * (((x+mu)^2 + y^2 + z^2)^(-3/2)   +...
                               -3*(x+mu)^2*(((x+mu)^2 + y^2 + z^2))^(-5/2)) +...
                   -    mu  * (((x+mu-1)^2 + y^2 + z^2)^(-3/2) +...
                               -3 * (x+mu-1)^2*(((x+mu-1)^2 + y^2 + z^2))^(-5/2));
dg1dy = @(x,y,z)    3*(1-mu) * y *   (x+mu) * ((x+mu)^2 + y^2 + z^2)^(-5/2) +...
                  + 3*    mu * y * (x+mu-1) * ((x+mu-1)^2 + y^2 + z^2)^(-5/2);
dg1dz = @(x,y,z)    3*(1-mu) * z *   (x+mu) * ((x+mu)^2 + y^2 + z^2)^(-5/2) +...
                  + 3*    mu * z * (x+mu-1) * ((x+mu-1)^2 + y^2 + z^2)^(-5/2);
%---------------------------------- g2 ------------------------------------
dg2dx = @(x,y,z) dg1dy(x,y,z);
dg2dy = @(x,y,z) 1 - (1-mu) * (((x+mu)^2 + y^2 + z^2)^(-3/2)   +...
                               -3*y^2*(((x+mu)^2 + y^2 + z^2))^(-5/2)) +...
                   -    mu  * (((x+mu-1)^2 + y^2 + z^2)^(-3/2) +...
                               -3 * y^2*(((x+mu-1)^2 + y^2 + z^2))^(-5/2));
dg2dz = @(x,y,z)    3*(1-mu) *        z * y * ((x+mu)^2 + y^2 + z^2)^(-5/2) +...
                  + 3*    mu *        z * y * ((x+mu-1)^2 + y^2 + z^2)^(-5/2);
%---------------------------------- g3 ------------------------------------
dg3dx = @(x,y,z) dg1dz(x,y,z);
dg3dy = @(x,y,z) dg2dz(x,y,z);
dg3dz = @(x,y,z) 0 - (1-mu) * (((x+mu)^2 + y^2 + z^2)^(-3/2)   +...
                               -3 * z^2*(((x+mu)^2 + y^2 + z^2))^(-5/2)) +...
                   -    mu  * (((x+mu-1)^2 + y^2 + z^2)^(-3/2) +...
                               -3 * z^2*(((x+mu-1)^2 + y^2 + z^2))^(-5/2));
%SubMatricies
G = @(x,y,z) [dg1dx(x,y,z), dg1dy(x,y,z), dg1dz(x,y,z);
              dg2dx(x,y,z), dg2dy(x,y,z), dg2dz(x,y,z);
              dg3dx(x,y,z), dg3dy(x,y,z), dg3dz(x,y,z)];
H = [0, 2, 0;
    -2, 0, 0;
     0, 0, 0];
%State Matrix
A = @(x,y,z) [O               I;
              G(x,y,z)        H];
%Make STM
STM  = reshape(s,6,6)';
dSTM = A(x(1),x(2),x(3)) * STM;
dstm = dSTM';
dstm = dstm(:);
%-------------------------------------------------------------------------
%The last 6 entries are the vector field for the 2-Body problem.
%These are stored in ’dx’.
r1 = sqrt((x(1)+mu)^2   + x(2)^2 + x(3)^2);
r2 = sqrt((x(1)+mu-1)^2 + x(2)^2 + x(3)^2);
%Set the derivatives of the state
dx = [x(4);
      x(5);
      x(6);
      x(1)+2*x(5)-((1-mu)*(x(1)+mu))/(r1^3)-mu*(x(1)-(1-mu))/(r2^3);
      x(2)-2*x(4)-(1-mu)/(r1^3)*x(2)-mu*(x(2))/(r2^3);
      -(1-mu)/(r1^3)*x(3)-mu*(x(3))/(r2^3)];
%-------------------------------------------------------------------------
% Put it all toghether and pass back to integrator
dy = [dx; dstm];
end