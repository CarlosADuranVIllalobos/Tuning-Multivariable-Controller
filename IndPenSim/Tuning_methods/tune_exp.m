function M = tune_exp(mina, maxa, u_length, num_u)
%Function to select the parameters of the diagonal suppresion matrix
%exponentially
%Copyright Carlos Alberto Duran Villalobos, February 2022-University of Manchester
     n = u_length;
     gf = (maxa/mina)^(1/(n - 1));
     k = 1:n;
     y = mina*gf.^(k-1);
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

