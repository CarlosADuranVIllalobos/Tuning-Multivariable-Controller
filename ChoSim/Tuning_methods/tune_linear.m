function M = tune_linear(mina, maxa, M_length, num_u)
%Function to select the parameters of the diagonal suppresion matrix linearly 
%Copyright Carlos Alberto Duran Villalobos, February 2022-University of Manchester
     n = M_length;
     gf = (maxa-mina)/n;
     k = 1:n;
     y = mina+gf*(k-1);
     M=ones(num_u*n,1);
     for i=1:num_u
         val=1;
         for k=1:n
             %index=(k-1)*num_x + (k-ustart)*num_u + i + num_x+invar;
             M(k+i*n-n)=y(val);
             val=val+1;
         end
     end
     
end

