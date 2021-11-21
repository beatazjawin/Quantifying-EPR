%Given two assemblages, check whether one can be converted to the other
%using local operations and shared randomness (LOSR)


function output=LOSR_conversions(d, Sig, SigT)

% Inputs:
% two assemblages - Sig and SigT
% Sig(d,d,A,X) is an assemblage with elements: Sig(:,:,a,x) = sigma_a|x (for |X|=|A|=2)
% d is the dimension of Bob's Hilbert space

% Output:
% Given two assemblages, it tells whether Sig can be LOSR-transformed to
% SigT (out=1) or not (out=0)

L=[];
for Q=1:2
for W=1:2
for E=1:2
for R=1:2
for T=1:2
for Y=1:2
    L=[L; Q W E R T Y];
end
end
end
end
end
end

 cvx_begin SDP quiet
    cvx_solver sdpt3 
    cvx_precision default
        variable W(d*d,d*d,64) complex % the Choi state W over subsystems B' and B
        expressions aux(d,d,2,2) complex 
        expressions auxn(d,d) complex 
        expressions ptr(d,d,64) complex

  % conditions on the hidden states

        for lambda=1:64
            W(:,:,lambda) == hermitian_semidefinite(d*d);
        end

        
        auxn=zeros(d,d);
        for lambda=1:64
            auxn=auxn+PartialTrace(W(:,:,lambda),1);
        end
        auxn==eye(d)/d; 


        for lambda=1:64
            ptr(:,:,lambda)=PartialTrace(W(:,:,lambda),1);
            for i=1:2
                for j=1:2
                    if (i~=j)
                        ptr(i,j,lambda)==0;
                    end
                end
            end
        ptr(1,1,lambda)==ptr(2,2,lambda);
                        
        end
        
      
       
   
  % check if Sig can be converted to SigT with LOSR operations
  
        for xp=1:2
            for ap=1:2
                 aux(:,:,ap,xp)=zeros(d,d);
                 for lambda=1:64
                       for x=1:2
                           for a=1:2
                                aux(:,:,ap,xp) = aux(:,:,ap,xp) + (ap==f(xp,a,L(lambda,:)))*(x==g(xp,L(lambda,:)))*d*PartialTrace(W(:,:,lambda)*Tensor(eye(d),Sig(:,:,a,x)'));
                           end
                       end
                 end
                 SigT(:,:,ap,xp) == aux(:,:,ap,xp);
            end
        end

               
    cvx_end
    
    output=1-min(1,cvx_optval);
    
end


function output=g(xp,L)
    output = mod(L(5)*xp + L(6), 2)+1 ;
end

function output=f(xp,a,L)
    output = mod(L(1)*a*xp+L(2)*a+L(3)*xp+L(4), 2)+1 ;
end 