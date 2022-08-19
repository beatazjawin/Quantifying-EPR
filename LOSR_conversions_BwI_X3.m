%Given two Bob-with-Input assemblages, check whether one can be converted to the other
%using local operations and shared randomness (LOSR)

%This code is for |X|=3, |A|=|Y|==2

function output=LOSR_conversions_BwI_X3(d, Sig, SigT)


% Inputs:
% d is the dimension of Bob's Hilbert space
% Sig(d,d,A,X,Y) is an assemblage with elements: Sig(:,:,a,x,y) = sigma_a|xy

% Output:
% Given two channel assemblages, it tells whether Sig can be LOSR-transformed to
% SigT (out=1) or not (out=0)

L=[];
for Q=1:3
for W=1:3
for E=1:3
for R=1:2
for T=1:2
for Y=1:2
for U=1:2
for O=1:2
for P=1:2
for S=1:2
for F=1:2
    L=[L; Q W E R T Y U O P S F];
end
end
end
end
end
end
end
end
end
end
end

 cvx_begin SDP quiet
 cvx_solver sedumi
 cvx_precision high    
        variable W(d*d,d*d,6912,2) complex
        expressions aux(d,d,2,3,2) complex 
        expressions auxn(d,d) complex 
        expressions ptr(d,d,6912,2) complex

% conditions on the map W
 

        for lambda=1:6912
            for yp=1:2
                W(:,:,lambda,yp) == hermitian_semidefinite(d*d);
            end
        end
        
        
        for yp=1:2
            auxn=zeros(d,d);
            for lambda=1:6912
                auxn=auxn+PartialTrace(W(:,:,lambda,yp),1);
            end
            auxn==eye(d)/d;
        end
         
        
        
        for lambda=1:6912
            for yp=1:2
                ptr(:,:,lambda,yp)=PartialTrace(W(:,:,lambda,yp),1);
                for i=1:2
                    for j=1:2
                        if (i~=j)
                            ptr(i,j,lambda,yp)==0;
                        end
                    end
                end
            ptr(1,1,lambda,yp)==ptr(2,2,lambda,yp);
            end            
        end
        
        for lambda=1:6912
            PartialTrace(W(:,:,lambda,1),1)==PartialTrace(W(:,:,lambda,2),1); 
        end
        
  % check if Sig can be converted to SigT with LOSR operations  
  % the deterministic strategies in L are labeled as L[]=[ x(x'=1), x(x'=2), x(x'=3), a'(a=1,x'=1), a'(a=1,x'=2), a'(a=1,x'=3),
  %a'(a=2,x'=1), a'(a=2,x'=2), a'(a=2,x'=3), y(y'=1), y(y'=2)]
       
        for xp=1:3
            for ap=1:2
                for yp=1:2
                    aux(:,:,ap,xp,yp)=zeros(d,d);
                    for lambda=1:6912
                           for x=1:3
                               for a=1:2
                                   for y=1:2
                                        aux(:,:,ap,xp,yp) = aux(:,:,ap,xp,yp) + (ap==L(lambda,a*3+xp))*(x==L(lambda,xp))*(y==L(lambda,12-yp))*d*PartialTrace(W(:,:,lambda,yp)*Tensor(eye(d),Sig(:,:,a,x,y).'));
                           
                                   end
                                end
                           end
                    end
 
                    SigT(:,:,ap,xp,yp) == aux(:,:,ap,xp,yp);
                end
            end
        end

               
    cvx_end
    
    output=1-min(1,cvx_optval);
    
end
