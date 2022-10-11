function [] = Main()
%use the modeset to pick either setup (1) for first run or something else
%for continuation of a previous run.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initiation and mode set
for dummy = 1:1;
clear all
mode = 1; %Mode 1 for setup, Mode 0 for continueing a previous run
tic
if mode == 1;
disp('We are started in setup mode!')
else
disp('We are started in continuation mode!')    
end
end

%Setup (create/load)
for dummy=1:1
if mode == 1
    NN = 18;
    beta = 5;
    %J=1
    tscram=0; %dummy
    NS=0; %dummy
    Z0=0; %dummy
    SN=6;
    ep=pi/8; %injection parameter
    g=pi/4; %coupling of double trace
    
    save('Setup','NN','beta','tscram','NS','SN','ep','g','Z0')
    
    %initialisation for experiments
    data1=[0 0];
    save('data1.mat','data1')
    
    data2=[0 0];
    save('data2.mat','data2')
    
    data3=[0 0];
    save('data3.mat','data3')
    
    data4=[0 0];
    save('data4.mat','data4')
    
    data5=[0 0];
    save('data5.mat','data5')
    
    data6=[0 0];
    save('data6.mat','data6')
    
else
    load('Setup.mat')
end
disp('setup done!')
t=toc;
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
end

%Fermions (Create) stored in sparse array X
for dummy=1:1
   %splitting of strange dimensions that cannot be constructed from the
%recursion.
if NN<1
    error('Dimensions do not make sense!')
elseif NN==1
    GA{1}=1;    
else

%generating the matrixes for NN=2, which is used for the basis of
%iterations.
n=floor(NN/2);
od=NN-2*n;

G{1}=[0,1;1,0];
G{2}=[0,-1i;1i,0];

GA=G;

%now iterate for the even dimensions.
for n1=2:n
    
   %generate 2*a1-2 matrices from previous step. but first generate the
   %zero blocks and identity blocks used for constructing the new matrcies.
   f0 = sparse(2^(n1-1),2^(n1-1));
   h0 = sparse(2^(n1-2),2^(n1-2));
   f1 = speye(2^(n1-1));
   h1 = speye(2^(n1-2));
   
   for a=1:2*n1-2
       G{a}=sparse([f0,GA{a};GA{a},f0]);    
   end
   
   %now generate the new matrixes.
   Gd=[h1,h0;h0,-h1];
   G{a+1}=[f0,Gd;Gd,f0];
   G{a+2}=[f0,-1i*f1;1i*f1,f0];
   
   GA=G;
end

%add one more matrix for odd dimensions.
if od==1
   f0 = sparse(2^(n1-1),2^(n1-1));
   f1 = speye(2^(n1-1));
   GA{NN}=[f1,f0;f0,-f1];
end

 

end

for h =1:NN
X{h}=GA{h}./sqrt(2);
end
clear f0 f1 G GA h0 h1 od n   
idd=speye(2^(floor(NN/2)));

disp('fermions done!')
t=toc;
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
end

%Hamiltonian (create/load), eigensystem (create)
for dummy=1:1
   if mode == 1
       H=zeros(2^(floor(NN/2)));   %starting point for the hamiltonian.

%adjust the hamiltonian for all the interactions.
for t1=1:(NN-3)
    for t2=(t1+1):(NN-2)
        for t3=(t2+1):(NN-1)
            for t4=(t3+1):NN
                J1=normrnd(0,sqrt(6/NN^3));
                H=H+J1*X{t1}*X{t2}*X{t3}*X{t4};
            end
        end
    end
end
    
    save('Hamiltonian.mat','H')
   else
       load('Hamiltonian.mat')
   end
   [V,E]=eig(H);  %V.D.inv(V) = H
   [~,permutation]=sort(diag(E));
   E=E(permutation,permutation);
   V=V(:,permutation); 
   
   
   
   Ediag = sparse(E);
   Evec = diag(E);  %colomn vector, smallest to biggest
   
   Z0 = trace(expm(-beta*Ediag));
   E0 = trace(expm(-beta*Ediag)*Ediag)/Z0;
   
   rc = expm(beta*Ediag/2);
   lc = expm(-beta*Ediag/2);
   
   disp('Hamiltonian done!')
t=toc;
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
end
    
%state (create/load)
for dummy=1:1
   if mode == 1
      
       
      n=length(Evec);
       
      Estate=Evec*0;
      n0=find(E0-Evec<0,1);
      Estate(n0)=randn + 1i*randn;
      
      NS = floor(2^(NN/2)/(Evec(n) - Evec(1))*0.015) + 3;
      

      
      
      upcounter = 0;
      downcounter = 0;
      
      
      for ns = 2:NS-1
          if dot(Estate,Ediag*Estate)/dot(Estate,Estate) > E0
              downcounter=downcounter+1;
              Estate(n0-downcounter)=randn+1i*randn;
          else
              upcounter=upcounter+1;
              Estate(n0+upcounter)=randn+1i*randn;
          end
      end
      
      if dot(Estate,Ediag*Estate)/dot(Estate,Estate) > E0
          ofset=-downcounter-1;
      else
          ofset=upcounter+1;
      end
      
      Estate(n0+ofset)= sqrt(abs((E0*dot(Estate,Estate)-real(dot(Estate,Ediag*Estate)))/(Evec(n0 + ofset) -E0)));
      Estate=Estate/norm(Estate);
      
           
      
      save('state.mat', 'Estate')
      
      F = -logm(Z0)./beta;
      S=beta*(E0-F);
      
      tscram= beta*log(S)/(2*pi);
          
      
      save('Setup.mat','NN','beta','tscram','NS','SN','ep','g','Z0') %appending tscram, and NS
       
   else
       
       load('state.mat')
       
       
   end
    

    disp('state done!')
t=toc;
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    
end
    
%going to energy state fermions XE= V.X.V^-1 
for dummy=1:1
%and simmilar for the spins S=2i X X
%Only doing this for fermion 1 and 2, and spins up to SN (setup)
    
iV=inv(V);
    for h =1:2
        XE{h}=V*X{h}*iV;
    end
    for h =1:SN
        SE{h}=2i*V*X{2*h-1}*X{2*h}*iV;
        SX{h}=2i*idd*X{2*h-1}*X{2*h};
    end
end
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Exp1: States in energy space, fermion and spin insertions
for dummy = 1:1
    load('data1.mat')
    if data1(1,1)==0
        

    
    
    EstateFint = lc*(cos(ep/sqrt(2))*idd+1i*sin(ep/sqrt(2))*XE{1})*rc*Estate;
    EstateFext = (cos(ep/sqrt(2))*idd+1i*sin(ep/sqrt(2))*XE{1})*Estate;
    EstateSint = lc*(cos(ep)*idd+1i*sin(ep)*SE{1})*rc*Estate;
    EstateSext = (cos(ep)*idd+1i*sin(ep)*SE{1})*Estate;    
    
    EFI = dot(EstateFint,Ediag*EstateFint);
    EFE = dot(EstateFext,Ediag*EstateFext);
    ESI = dot(EstateSint,Ediag*EstateSint);
    ESE = dot(EstateSext,Ediag*EstateSext);
    
    data1 = [Evec,Estate,EstateFint,EstateFext,EstateSint,EstateSext];
    save('data1.mat','data1')
    csvwrite('exp1R.csv',real(data1))
    csvwrite('exp1I.csv',imag(data1))
    
    exp1E =real([E0, EFI, EFE, ESI, ESE]);
    
    save('data1energy.mat','exp1E')
    csvwrite('exp1energy.csv',exp1E)
    
    disp('exp1 done!')
    t=toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    
    clear EstateFint EstateFext EstateSint EstateSext
    else
    disp('exp1 was already done!')
    t=toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    end

end

%Exp2: Matrix elements of a fermion probe around scrambling time
for dummy = 1:1
    load('data2.mat')
    if data2(1,1)==0
        
   n=2^floor(NN/2);     
        
   help1=spdiags(exp(1i*Evec*tscram),0,n,n)*XE{1}*spdiags(exp(-1i*Evec*tscram),0,n,n);
   probe = (help1*XE{2}+XE{2}*help1)^2;
   n=length(probe);
   Eprobe = diag(probe);
   
   help1=probe;
   help1(:,1)=[];
   help1(n,:)=[];
   Eup = [diag(help1);0];
   
   probe(1,:)=[];
   probe(:,n)=[];
   Edo =  [diag(help1);0];
   
   can=exp(-beta*Evec);
   
   data2=full([Evec,Estate,Eprobe,Eup,Edo,can]); 
   
   save('data2.mat','data2')
   csvwrite('exp2R.csv',real(data2))
   csvwrite('exp2I.csv',imag(data2))
    

    
    disp('exp2 done!')
    t=toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    clear probe Uup Edo can help1
    else
    disp('exp2 was already done!')
    t=toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    end
end

%Exp3: One point functions in different back grounds    
for dummy = 1:1
    tstart = -4*tscram;
    tend = 4*tscram;
    dstep = tscram/25;
    steps = (tend-tstart)/dstep;
    load('data3.mat')
    
    if data3(end,1) == steps  
        disp('exp3 was already done!')
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS')) 
    else
    start = data3(end,1)+1;
    
    n=2^floor(NN/2);
    
    EstateFext = (cos(ep/sqrt(2))*idd+1i*sin(ep/sqrt(2))*XE{1})*Estate;
    EstateSext = (cos(ep/sqrt(2))*idd+1i*sin(ep/sqrt(2))*SE{1})*Estate; 
    
    for step = start:steps
        
    data3(step,1)=step;    
    ts = step*dstep+tstart;
    data3(step,2)=ts;
    
    help=spdiags(exp(1i*Evec*ts),0,n,n)*XE{1}*spdiags(exp(-1i*Evec*tscram),0,n,n);
    data3(step,3) = dot(Estate,help*Estate);
    data3(step,4) = dot(EstateFext,help*EstateFext);
        
    help=spdiags(exp(1i*Evec*ts),0,n,n)*SE{1}*spdiags(exp(-1i*Evec*ts),0,n,n);
    data3(step,5) = dot(Estate,help*Estate);
    data3(step,6) = dot(EstateSext,help*EstateSext);   
        
    save('data3.mat','data3')        
        
    end
    
    csvwrite('exp3R.csv',real(data3))
    csvwrite('exp3I.csv',imag(data3))
    
disp('exp3 done!')
t=toc;
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
clear EstateFext EstateSext help
    end
end
    
%Exp4: (anti-)commutator in canonical and state  
for dummy = 1:1
    tstart = 0;
    tend = 4*tscram;
    dstep = tscram/40;
    steps = (tend-tstart)/dstep;
    load('data4.mat')
    
    if data4(end,1) == steps  
        disp('exp4 was already done!')
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS')) 
    else
       start = data4(end,1)+1;
        
       n=2^floor(NN/2);

    for step = start:steps
        
    data4(step,1)=step;    
    ts = step*dstep+tstart;
    data4(step,2)=ts;
    boltz = expm(-beta*Ediag)/Z0;
    
    
    help1=spdiags(exp(1i*Evec*ts),0,n,n)*XE{1}*spdiags(exp(-1i*Evec*ts),0,n,n);
    help = (help1*XE{2}+XE{2}*help1)^2;
    data4(step,3) = dot(Estate,help*Estate);
    data4(step,4) = trace(boltz*help);
        
    help1=spdiags(exp(1i*Evec*ts),0,n,n)*SE{1}*spdiags(exp(-1i*Evec*ts),0,n,n);
    help = -(help1*SE{2}-SE{2}*help1)^2;
    data4(step,5) = dot(Estate,help*Estate);
    data4(step,6) = trace(boltz*help);   
        
    save('data4.mat','data4')        
        
    end
    
    csvwrite('exp4R.csv',real(data4))
    csvwrite('exp4I.csv',imag(data4))
    
    disp('exp4 done!')
    t=toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    clear help help1 
    end
    
end
  
%Exp5: two-point function (with and without analyticly continued second operator)
for dummy = 1:1
    tstart = -4*tscram;
    tend = 4*tscram;
    dstep = tscram/25;
    steps = (tend-tstart)/dstep;
    load('data5.mat')
    
    if data5(end,1) == steps  
        disp('exp5 was already done!')
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS')) 
    else
        start = data5(end,1)+1;
        
        n=2^floor(NN/2);
    
    help1 = lc*XE{1}*rc;
    help2 = lc*SE{1}*rc;
    boltz = expm(-beta*Ediag)/Z0;
    
    for step = start:steps
        
    data5(step,1)=step;    
    ts = step*dstep+tstart;
    data5(step,2)=ts;
    
    
    
    helpt=spdiags(exp(1i*Evec*ts),0,n,n)*XE{1}*spdiags(exp(-1i*Evec*ts),0,n,n);
    help = helpt*XE{1};
    data5(step,3) = dot(Estate,help*Estate);
    data5(step,4) = trace(boltz*help);
    help = helpt*help1;
    data5(step,5) = dot(Estate,help*Estate);
    data5(step,6) = trace(boltz*help);
        
    help1=spdiags(exp(1i*Evec*ts),0,n,n)*SE{1}*spdiags(exp(-1i*Evec*ts),0,n,n);
    help = help1*SE{1};
    data5(step,7) = dot(Estate,help*Estate);
    data5(step,8) = trace(boltz*help);   
    help = help1*help2;
    data5(step,9) = dot(Estate,help*Estate);
    data5(step,10) = trace(boltz*help);
        
    save('data5.mat','data5')        
        
    end
    
    csvwrite('exp5R.csv',real(data5))
    csvwrite('exp5I.csv',imag(data5))
    
disp('exp5 done!')
t=toc;
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
clear help help1 helpt help2
    end
end    
    
%Exp6: double trace experiment
for dummy = 1:1
    
    
    
    
    tstart = -2*tscram;
    tend = 3*tscram;
    dstep = tscram/30;
    steps = (tend-tstart)/dstep;
    load('data6.mat')
    
    if data6(end,1) == steps  
        disp('exp6 was already done!')
        t=toc;
        disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS')) 
    else
        start = data6(end,1)+1;
        
    n=2^floor(NN/2);
    boltz = expm(-beta*Ediag)/Z0;   
    
    help1=spdiags(exp(1i*Evec*tscram),0,n,n)*SE{1}*spdiags(exp(-1i*Evec*tscram),0,n,n);
    EstateMod = lc;
    XstateMod = iV*EstateMod*(cos(ep)*idd+1i*sin(ep)*help1)*V;
    help2 = iV*lc*V;
    
    
    K=SN-1;
    
    for h = 2:SN
        XstateMod = cos(g/K)*XstateMod+1i*sin(g/K)*SX{h}*XstateMod*SX{h};
        help2 = cos(g/K)*help2+1i*sin(g/K)*SX{h}*help2*SX{h};
    end
     
    XSM1 = XstateMod;
    XSM2 = help2;
    
    
    for step = start:steps
        
    XstateMod = XSM1;
    help2 = XSM2;  
        
    data6(step,1)=step;    
    ts = step*dstep+tstart;
    data6(step,2)=ts;

    help3=iV*spdiags(exp(1i*Evec*ts),0,n,n)*SE{1}*spdiags(exp(-1i*Evec*ts),0,n,n)*V;
        
    XstateMod = help3*XstateMod;
    
    for h = 2:SN
        XstateMod = cos(g/K)*XstateMod-1i*sin(g/K)*SX{h}*XstateMod*SX{h};
    end    
    
    EstateMod = V*XstateMod*iV*(cos(ep)*idd-1i*sin(ep)*help1)*rc;
    
    data6(step,3) = dot(Estate,EstateMod*Estate);
    data6(step,4) = trace(boltz*EstateMod);  
    
    %and now without probe
        
    XstateMod = help3*help2;
    
    for h = 2:SN
        XstateMod = cos(g/K)*XstateMod-1i*sin(g/K)*SX{h}*XstateMod*SX{h};
    end    
    
    EstateMod = V*XstateMod*iV*rc;
    
    data6(step,5) = dot(Estate,EstateMod*Estate);
    data6(step,6) = trace(boltz*EstateMod); 
        
    save('data6.mat','data6')        
        
    end
    
    save('dpertstate.mat','XSM2')
    
    csvwrite('dpertstateR.csv',real(XSM2))
    csvwrite('dpertstateI.csv',imag(XSM2))
    
    csvwrite('exp6R.csv',real(data6))
    csvwrite('exp6I.csv',imag(data6))
    
    disp('exp6 done!')
    t=toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    end
end     
    
    
end


