% Problema parabolico: du/dt = c*d^2u/dx^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function heatEquationSolve(n,m)
  c = 0.25;
  % intervalos x y t
  x0 = 0; xL = 2; t0 = 0; tf = 5.5;
  % malla
  h = (xL-x0)/(n+1); 
  k = (tf-t0)/m;
  % Condiciones iniciales
  
  for i = 1:n;
    x(i) = i*h; 
    U(i) = -x(i)*(x(i)-2);%sin(pi*x(i));
  end
  % Valores en la frontera
  for i = 1:m+1;
    t(i)= t0 + i*k;
    f0(i)=0;     
    fL(i)=0;
  end
  MU = [ f0(1) U fL(1) ];  MUex = MU;
  A =zeros(n); 
  b=zeros(n,1);
  for i = 1:n-1;
    A(i,i)= -2; 
    A(i,i+1) = 1; 
    A(i+1,i) = 1;
  end
  A(n,n) = -2;
  A 
  r = c*k/h^2;
  I = eye(n);
  for i=1:m
    b(1) = (f0(i+1) + f0(i))/2;
    b(n) = (fL(i+1) + fL(i))/2;
    U = (inv(I-r*A/2)*((I+r*A/2)*U' + r*b ))';
    % MU contiene las soluciones aproximada por filas
    MU = [ MU ; f0(i+1) U fL(i+1) ];
    % Solucion exacta
    for j=1:n
      uex(j)= F(x(j),t(i));
    end
    MUex = [MUex ;f0(i+1) uex fL(i+1) ];
  end
  MU
  figure(1);
  surf(MU)
  title({'Ecuacion de Calor : Metodo Clanck Nicholson';['{Error maximo} = ',num2str(max(max(abs(MU-MUex))))]})
  %subplot (2, 1, 2)
  figure(2);
  surf(MUex)
  title({'Ecuacion de Calor F(x,t) = exp(-(pi^2)*t/4)*sin(pi*x)'})
end

function [z] = F(x,t)
  z = exp(-(pi^2)*t/4)*sin(pi*x);
end