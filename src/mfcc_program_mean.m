%Name:Muhammad Zubair SC
%Title: Automatic Sound Reconigtion System

function m = mfcc_program_mean(audp, fs)

N=length(audp);
preEmph=2* [1 -1];
aud = filter(preEmph,1,audp);
frameduration=0.025;
framelen=round(frameduration*fs);
overlap=round(0.01*fs);
lasframe=floor(N/framelen);
audf=zeros(framelen,lasframe);
filterd=zeros(framelen,lasframe);
nfft=nextpow2(framelen);
XFTful=zeros(2^nfft,lasframe);
fftlen=2^nfft/2;
Xft=zeros(fftlen,lasframe);
F=26;
M=zeros(F,fftlen);
k=0:fftlen-1;
wn=hamming(framelen);

 for o=1:lasframe
 audf(:,o)=aud(((o-1)*(framelen-overlap)+1):(o*framelen-(o-1)*overlap));
 filterd(:,o)=wn.*audf(:,o);
 XFTful(:,o)=abs(fft(filterd(:,o),2^nfft));
 Xft(:,o)=XFTful(1:(2^nfft)/2,o);
 end
 
 %% Filter Design
 lfk=k*fs/framelen;
 lfmax=fs;
 lfmin=100;
 phimax=2595*log10(1+lfmax/700);
 phimin=2595*log10(1+lfmin/700);
 Dphif=(phimax-phimin)/(F+1);
 phifC=(1:F)*Dphif;
 lfc=700*(10.^(phifC/2595)-1);

 for m=2:F-1;
     for k=1:(2^nfft)/2
 if((lfk(k)<lfc(m-1))||(lfk(k)>=lfc(m+1)))
     M(m,k)=0;
elseif ((lfc(m-1)<=lfk(k))&&(lfk(k)<lfc(m)))
M(m,k)=(lfk(k)-lfc(m-1))/(lfc(m)-lfc(m-1));
elseif ((lfc(m)<=lfk(k))&&(lfk(k)<lfc(m+1)))
M(m,k)=(lfk(k)-lfc(m+1))/(lfc(m)-lfc(m+1));;
end
    end
 end
LogEnergy=zeros(F,lasframe);
Ener=zeros((2^nfft)/2,1);
Energy=zeros(F,lasframe);

for o=1:lasframe
    for m=1:F;
    for n=1:(2^nfft)/2
 Ener(n)=M(m,n)*Xft(n,o);
    end
    Energy(m,o)=sum(Ener);
    Ener=zeros((2^nfft)/2,1);
    end
end
LogEnergy=log(Energy);

for o=1:lasframe
    for m=1:26
     if((LogEnergy(m,o)==Inf)||(LogEnergy(m,o)==-Inf)||isnan(LogEnergy(m,o)))
         LogEnergy(m,o)=0;
     end
    end
end
cn=dct(LogEnergy);

cnoverall=sum(cn');
cnoverall=cnoverall(1,1:13)';
m=cnoverall/lasframe;
