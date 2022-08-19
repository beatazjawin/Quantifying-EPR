%Given two Bob-with-Input assemblages, check whether one can be converted to the other
%using local operations and shared randomness (LOSR)

%This code is for |X|=|A|=|Y|=2

function output=LOSR_conversions_BwI_X2(d, Sig, SigT)

% Inputs:
% d is the dimension of Bob's Hilbert space
% Sig(d,d,A,X,Y) is an assemblage with elements: Sig(:,:,a,x,y) = sigma_a|xy

% Output:
% Given two channel assemblages, it tells whether Sig can be LOSR-transformed to
% SigT (out=1) or not (out=0)

L=[];
for Q=1:2
for W=1:2
for E=1:2
for R=1:2
for T=1:2
for Y=1:2
for U=1:2
for O=1:2
    L=[L; Q W E R T Y U O];
end
end
end
end
end
end
end
end

    cvx_begin SDP quiet
    cvx_solver SeDuMi
    cvx_precision best    
        variable W(d*d,d*d,256,2) complex 
        expressions aux(d,d,2,2,2) complex 
        expressions auxn(d,d) complex
        expressions ptr(d,d,256,2) complex

  % conditions on the map W

        for lambda=1:256
            for yp=1:2
                W(:,:,lambda,yp) == hermitian_semidefinite(d*d);
            end
        end
        
        for yp=1:2
            auxn=zeros(d,d);
            for lambda=1:256
                auxn=auxn+PartialTrace(W(:,:,lambda,yp),1);
            end
            auxn==eye(d)/d;
        end
         
        for lambda=1:256
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
        
        for lambda=1:256
            PartialTrace(W(:,:,lambda,1),1)==PartialTrace(W(:,:,lambda,2),1); 
        end
        
  % check if Sig can be converted to SigT with LOSR operations  
       
        for xp=1:2
            for ap=1:2
                for yp=1:2
                    aux(:,:,ap,xp,yp)=zeros(d,d);
                    for lambda=1:256
                           for x=1:2
                               for a=1:2
                                   for y=1:2
                                        aux(:,:,ap,xp,yp) = aux(:,:,ap,xp,yp) + (ap==f(xp,a,L(lambda,:)))*(x==g(xp,L(lambda,:)))*(y==h(yp,L(lambda,:)))*d*PartialTrace(W(:,:,lambda,yp)*Tensor(eye(d),Sig(:,:,a,x,y).'));
                           
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


function output=h(xp,L)
    output = mod(L(7)*xp + L(8), 2)+1 ;
end

function output=g(xp,L)
    output = mod(L(5)*xp + L(6), 2)+1 ;
end

function output=f(xp,a,L)
    output = mod(L(1)*a*xp+L(2)*a+L(3)*xp+L(4), 2)+1 ;
end 