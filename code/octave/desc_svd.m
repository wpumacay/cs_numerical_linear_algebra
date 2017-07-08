A=[1 1;
   1 0;
   0 1];
   
t=0:0.01:2*pi;

x=cos(t);
y=sin(t);
X=[x;y];XX=[];

subplot(1,2,1)
plot(x,y),grid,axis square

for i=1:length(x)
 x1=X(:,i);
 Xi=A*x1;
 XXX=[XX Xi];
 XX=XXX;
end

subplot(1,2,2)
plot3(XX(1,:),XX(2,:),XX(3,:),'r'),grid, axis square

AA=A'*A
[V,D]=eig(AA);
v1=V(:,1);v2=V(:,2);
Av1=A*v1, Av2=A*v2
norm(Av1),norm(Av2)