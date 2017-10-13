%Gram-Schmidt modificado
function [q,r]=gramsm(a)
[m n]=size(a);
for j=1:n
r(j,j)=norm(a(:,j));
a(:,j)=a(:,j)/r(j,j);
r(j,j+1:n)=a(:,j)'*a(:,j+1:n);
a(:,j+1:n)=a(:,j+1:n)-a(:,j)*r(j,j+1:n);
end
q=a;