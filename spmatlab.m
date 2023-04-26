pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

Nx=50;
Ny=50;
J=1;
W=1;
K=100;
n=5000;
tic;
id=zeros(1,Nx*Ny);
for i=1: Nx*Ny
    id(i)=1;
end
Id=diag(id);

H1c=zeros(Nx *Ny);
for i=0:Nx*Ny-1
    for j=0:Nx*Ny-1
        if (Ny==mod(i-j,Nx*Ny)) && (mod((j-mod(j,Ny))/Ny+mod(j,Ny),2)==0)
            H1c(i+1,j+1)=-J;
        elseif (Ny==mod(j-i,Nx*Ny)) && (mod((i-mod(i,Ny))/Ny+mod(i,Ny),2)==0)
            H1c(i+1,j+1)=-J;
        else
            H1c(i+1,j+1)=0;
        end
    end
end
H2c=zeros(Nx *Ny);
for i=0:Nx*Ny-1
    for j=0:Nx*Ny-1
        if (1==mod(i-j,Ny)) && (mod((j-mod(j,Ny))/Ny+mod(j,Ny),2)==0) && ((j-mod(j,Ny))/Ny==(i-mod(i,Ny))/Ny) && (mod(j,Ny)~= Ny-1)
            H2c(i+1,j+1)=-J;
        elseif (1==mod(j-i,Ny)) && (mod((i-mod(i,Ny))/Ny+mod(i,Ny),2)==0) && ((j-mod(j,Ny))/Ny==(i-mod(i,Ny))/Ny) && (mod(j,Ny)~= 0)
            H2c(i+1,j+1)=-J;
        else
            H2c(i+1,j+1)=0;
        end
    end
end
H3c=zeros(Nx *Ny);
for i=0:Nx*Ny-1
    for j=0:Nx*Ny-1
        if (Ny==mod(j-i,Nx*Ny)) && (mod((j-mod(j,Ny))/Ny+mod(j,Ny),2)==0)
            H3c(i+1,j+1)=-J;
        elseif (Ny==mod(i-j,Nx*Ny)) && (mod((i-mod(i,Ny))/Ny+mod(i,Ny),2)==0)
            H3c(i+1,j+1)=-J;
        else
            H3c(i+1,j+1)=0;
        end
    end
end
H4c=zeros(Nx *Ny);
for i=0:Nx*Ny-1
    for j=0:Nx*Ny-1
        if (1==mod(j-i,Ny)) && (mod((j-mod(j,Ny))/Ny+mod(j,Ny),2)==0) && ((j-mod(j,Ny))/Ny==(i-mod(i,Ny))/Ny) && (mod(j,Ny)~= 0)
            H4c(i+1,j+1)=-J;
        elseif (1==mod(i-j,Ny)) && (mod((i-mod(i,Ny))/Ny+mod(i,Ny),2)==0) && ((j-mod(j,Ny))/Ny==(i-mod(i,Ny))/Ny) && (mod(j,Ny)~= Ny-1)
            H4c(i+1,j+1)=-J;
        else
            H4c(i+1,j+1)=0;
        end
    end
end

loc0=zeros(K,Nx/2,Ny/2);

parfor k=1:K
disv = zeros(1,Nx*Ny);
for i = 1:Nx*Ny
       disv(i) = W *(2*rand-1);
end 
dis=diag(disv);


H1=H1c+dis;
H2=H2c+dis;
H3=H3c+dis;
H4=H4c+dis;
H5=dis;

M=expm((1i*pi/2).*H5)*expm((1i*pi/2).*H4)*expm((1i*pi/2).*H3)*expm((1i*pi/2).*H2)*expm((1i*pi/2).*H1);
Mn=M^n;

loc1=zeros(Nx/2,Ny/2);
for i=0:Nx/2-1
    for j=0:Ny/2-1
        loc1(i+1,j+1)=abs(Mn(Ny *mod(Nx/2+i,Nx)+mod(Ny/2+j,Ny),Ny* mod(Nx/2,Nx)+mod(Ny/2,Ny))/Mn(Ny* mod(Nx/2,Nx)+mod(Ny/2,Ny),Ny* mod(Nx/2,Nx)+mod(Ny/2,Ny)));
    end
end
k

loc0(k,:,:)=loc1;
%loc0;


end
%loc0

loc=zeros(Nx/2,Ny/2);
for i=0:Nx/2-1
    for j=0:Ny/2-1
        vecloc=loc0(:,i+1,j+1);
        
loc(i+1,j+1)=max(vecloc);
    end
end
loc
csvwrite('localize.csv',loc)