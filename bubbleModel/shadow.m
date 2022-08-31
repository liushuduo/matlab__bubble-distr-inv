function flag=shadow(x1,R1,x2,R2)

if (R1<R2), flag=0; return, end

flag=norm(x1(:)-x2(:))<sqrt(R1);
   
