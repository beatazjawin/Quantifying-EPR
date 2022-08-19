%Given two channel assemblages, check whether one can be converted to the other
%using local operations and shared randomness (LOSR)

%This code is for |X|=3, |A|=2

function output=LOSR_conversions_channel_X3(din,dout,J1,J2)

% Inputs:
% din is the dimension of Bob's input Hilbert space
% dout is the dimension of Bob's output Hilbert space
% two assemblages in the Choi form - J1 and J2
% J(din*dout,din*dout,A,X) is an assemblage with elements: J(:,:,a,x) = J_a|x 

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
end
end
    cvx_begin SDP quiet
    cvx_solver sedumi
    cvx_precision high
        variable J(din*din*dout*dout,din*din*dout*dout,1728) complex % the Choi state of the LOSR map defined over subsystems C'BCB'
        variable F(din*din,din*din,1728) complex % the Choi state defined over subsystems CC'
        expressions auxn1 complex
        expressions auxn2 complex
        expressions ptr1(din*din,din*din,1728) complex
        expressions ptr2(din,din,1728) complex
        expressions aux(din*dout,din*dout,2,3) complex
      
 % conditions on the maps J and F
 
        for lambda=1:1728
           J(:,:,lambda) == hermitian_semidefinite(din*din*dout*dout);
        end
        
        for lambda=1:1728
            F(:,:,lambda) == hermitian_semidefinite(din*din);
        end
        
        auxn1=zeros(din*dout,din*dout);
        for lambda=1:1728
            auxn1=auxn1+PartialTrace(J(:,:,lambda),[3,4],[din,dout,din,dout]);
        end
        auxn1==eye(din*dout)/(din*dout); 

        auxn2=zeros(din,din);
        for lambda=1:1728
            auxn2=auxn2+PartialTrace(F(:,:,lambda),1,[din,din]);
        end
        auxn2==eye(din)/din; 

        
        for lambda=1:1728
            ptr1(:,:,lambda)=PartialTrace(J(:,:,lambda),[3,4],[din,dout,din,dout]);
            for i=1:din*dout
                for j=1:din*dout
                    if (i~=j)
                        ptr1(i,j,lambda)==0;
                    end
                end
                ptr1(1,1,lambda)==ptr1(i,i,lambda);
            end         
        end
        

        for lambda=1:1728
            ptr2(:,:,lambda)=PartialTrace(F(:,:,lambda),1,[din,din]);
            for i=1:2
                for j=1:2
                    if (i~=j)
                        ptr2(i,j,lambda)==0;
                    end
                end
                ptr2(1,1,lambda)==ptr2(i,i,lambda);
            end         
        end


        for lambda=1:1728
            ns=PermuteSystems(PartialTrace(J(:,:,lambda),4,[din,dout,din,dout]),[3,1,2],[din,dout,din]); 
            ns == Tensor(F(:,:,lambda),eye(dout))/dout;
        end  
        
  % check if J1 can be converted to J2 with LOSR operations  
  % the deterministic strategies in L are labeled as L[]=[ x(x'=1), x(x'=2), x(x'=3), a'(a=1,x'=1), a'(a=1,x'=2), a'(a=1,x'=3),
  %a'(a=2,x'=1), a'(a=2,x'=2), a'(a=2,x'=3)]
  
        for xp=1:3
            for ap=1:2
                 aux(:,:,ap,xp)=zeros(din*dout,din*dout);
                 for lambda=1:1728
                       for x=1:3
                           for a=1:2
                                insertensor=PermuteSystems(Tensor(eye(din*dout),PartialTranspose(J1(:,:,a,x),1,[din,dout])),[1,3,4,2],[din,dout,din,dout]); 
                                insert=insertensor*PartialTranspose(J(:,:,lambda),3,[din,dout,din,dout]); 
                                aux(:,:,ap,xp) = aux(:,:,ap,xp) + (ap==L(lambda,a*3+xp))*(x==L(lambda,xp))*din*dout*PermuteSystems(PartialTrace(insert,[2,3],[din,dout,din,dout]),[2,1],[din,dout]);
                           
                           end
                       end
                 end
                 J2(:,:,ap,xp) == aux(:,:,ap,xp);
            end
        end          
    cvx_end
   
    output=1-min(1,cvx_optval);

end
