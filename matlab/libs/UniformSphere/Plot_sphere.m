% Demonstrate the mySphere function
N = 1000; 
[X,Y,Z,N_new] = mySphere(N);

switch 2
    case 1
        % double the radius
        X = X.*2;
        Y = Y.*2;
        Z = Z.*2;
    case 2
        % offset sphere
        X = X + 1;
        Y = Y + 1;
        Z = Z + 1;        
end

figure
for i = 1:N_new
   hold on
   plot3(X(i),Y(i),Z(i),'b.')
end