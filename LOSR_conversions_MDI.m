%Given two measurement-device-independent assemblages, check whether one can be converted to the other
%using local operations and shared randomness (LOSR)

%This code is for |X|=3, |A|=|Y|==2

function output=MDI_conversions(d,J1,J2)

% Inputs:
% d is the dimension of Bob's Hilbert space
% two assemblages in the Choi form - J1 and J2
% J(d,d,A,B,X) is an assemblage with elements: J(:,:,a,b,x) = J_ab|x 

% Output:
% Given two channel assemblages, it tells whether J1 can be LOSR-transformed to
% J2 (out=1) or not (out=0)

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
    L=[L; Q W E R T Y U O P];
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
    cvx_precision best
        variable J(d*d,d*d,2,2,1728) complex % the Choi state of the LOSR map defined over subsystems In In'
        expressions auxn1 complex
        expressions auxn2 complex
        expressions auxn3 complex
        expressions aux(d,d,2,2,3) complex

  % conditions on the maps J      
        for b=1:2
            for lambda=1:1728
                for bp=1:2
                    J(:,:,bp,b,lambda) == hermitian_semidefinite(d*d);
                end
            end
        end
        
        for b=1:2
            auxn1=zeros(d,d);
            for bp=1:2
                for lambda=1:1728
                    auxn1=auxn1+PartialTrace(J(:,:,bp,b,lambda),1,[d,d]);
                end
            end
            auxn1==eye(d)/d;
        end


        for b=1:2
            for lambda=1:1728
                auxn2=zeros(d,d);
                for bp=1:2
                    auxn2=auxn2+PartialTrace(J(:,:,bp,b,lambda),1,[d,d]);
                end
                for i=1:d
                    for j=1:d
                        if (i~=j)
                            auxn2(i,j)==0;
                        end
                    end
                end
                auxn2(1,1)==auxn2(2,2); 
            end         
        end


        for b=1:2
            for lambda=1:1728
                auxn3=zeros(d*d,d*d);
                for bp=1:2
                    auxn3=auxn3+J(:,:,bp,b,lambda);
                end
                auxn3 == hermitian_semidefinite(d*d);
            end
        end
  
  % check if J1 can be converted to J2 with LOSR operations  
  % the deterministic strategies in L are labeled as L[]=[ x(x'=1), x(x'=2), x(x'=3), a'(a=1,x'=1), a'(a=1,x'=2), a'(a=1,x'=3),
  %a'(a=2,x'=1), a'(a=2,x'=2), a'(a=2,x'=3)]   
       
        for xp=1:3
            for ap=1:2
                for bp=1:2
                    aux(:,:,ap,bp,xp)=zeros(d);
                    for lambda=1:1728
                        for x=1:3
                            for a=1:2
                                for b=1:2
                                    link=d*PartialTrace(Tensor(J1(:,:,a,b,x),eye(2))*PartialTranspose(J(:,:,bp,b,lambda),1,[2,2]),1,[2,2]);
                                    aux(:,:,ap,bp,xp) = aux(:,:,ap,bp,xp) + (ap==L(lambda,a*3+xp))*(x==L(lambda,xp))*link;
                                end
                            end
                        end
                    end
                    J2(:,:,ap,bp,xp) == aux(:,:,ap,bp,xp);
                end
            end
        end 
        
      
    cvx_end
    
    output=1-min(1,cvx_optval);

    
end
