%% Creazione matrice e calcolo autovalori

clear all
close all 
clc

L=1.2; %[m]
h=8e-3; %[m]
b=0.04; %[m]
rho=2700; %[Kg/m^3]
E=68e9; %[Pa]

%moment of inertia
J= (b*h^(3))/12; %[m^4]
A=b*h; %[m^2]
m= rho*A; %[Kg/m]
%g=gamma
%w=omega
g=sym('g');

%%
%creazione matrice

H=[1 0 1 0; 0 1 0 1; -cos(g*L) -sin(g*L) cosh(g*L) sinh(g*L); sin(g*L) -cos(g*L) sinh(g*L) cosh(g*L)];

%studio determinante matrice
determinante= det(H);

%determinante numerico
%siamo interessati al range di frequenze 0-200Hz quindi omega varia tra
%0:1256,6
fmax=200;
f=[0 fmax];
w_in=2*pi.*f; %valori min e max di omega, il range di omega di cui siamo interessati
w_values=w_in(1):0.5:w_in(2);

fun=matlabFunction(determinante);% trasformo la funzione simbolica, così il codice è più veloce

det_num = zeros(length(w_values),1);

for i = 1:length(w_values)

    g=(m/(E*J)*w_values(i)^(2))^(1/4);
    det_num(i)=fun(g);
end

% Graph of the determinant as a function of omega 
figure(1)
plot(w_values, det_num, 'b',LineWidth=2)
grid on
xlabel('\omega: Angular speed')
ylabel('values of determinant of [H]')
hold on
plot(w_values,zeros(length(w_values),1),'k')
ylim([-2e4,2e4])

%zeros of the determinat that is a non linear function
%we want to compute det(H)=0, where det(H) is a non linear function. Quindi
%trovare i valori di omega per i quali la matrice H è singolare. Questi
%valori corrispondono alle frequenze naturali. Prima troviamo i valori di
%gamma in cui si annulla il determinante e da questi calcoliamo i valori di
%omega.

g=@(w) (m/(E*J)*w^(2))^(1/4);
g_min=g(w_in(1));
g_max=g(w_in(2));% g va da 0 a 10.4127
w=@(g) g^(2)*sqrt(E*J/m);
fun=matlabFunction(determinante);% da symbol diventa una funzione in funzione di g

%risolviamo la funzione col comando fsolve, nell'intervallo di interesse
%delle gamma.Essendo una funzione non lineare fsolve usa un metodo iterativo, quindi ha bisogno 
% di valori iniziali dai quali poter iniziare ad iterare.
%s1=round(fsolve(fun,[1:(g_max/2)]),4,'decimals');
%fsolve troverà più soluzioni, prendendo nelle diverse iterazioni valori
%iniziali diversi, compresi tra 1 e (g_max/2). Incrementando il valore iniziale di uno 
%ogni qual volta inizia un nuovo calcolo di uno zero.
% Non consideriamo il caso in cui gamma è zero perchè la
%soluzione è gamma uguale a zero e quindi omega uguale a zero. In questo
%caso il sistema è fermo e non ci interessa studiarlo. Per questo partiamo
%dal valore iniziale  1.
%dividiamo in due intervalli perchè altrimenti eccediamo il numero max di
%iterazioni

s1=round(fsolve(fun,1:(g_max/2)),4,'decimals'); %prendo le cifre significative fino alla 4 decimale dopo la virgola
s2=round(fsolve(fun,(g_max/2):g_max),4,'decimals');
sol_g= unique([s1 s2]);

% troviamo tutti i valori di gamma per i quali la matrice H è singolare 
% da qua ricaviamo i valori di omega, cioè le frequenze proprie.
w1=w(sol_g(1));
w2=w(sol_g(2));
w3=w(sol_g(3));
w4=w(sol_g(4));
w_nf=[w1 w2 w3 w4]; % [rad/s]
nat_f=w_nf./(2*pi); %frequenze naturali [Hz]

figure(2)
plot(w_values, det_num, 'b',LineWidth=2)
grid on
xlabel('\omega: Angular speed')
ylabel('values of determinant of [H]')
hold on
plot(w_values,zeros(length(w_values),1),'--k')
hold on
plot(w1,g(w1),'.r',MarkerSize=15)
hold on
plot(w2,g(w2),'.r',MarkerSize=15)
hold on
plot(w3,g(w3),'.r',MarkerSize=15)
hold on
plot(w4,g(w4),'.r',MarkerSize=15)
ylim([-2e4,2e4])
title('Graph of det[H]')
legend('det[H]','x axis','zeros')

%% Calcolo dei mode shapes
z=struct('z_m1',[],'z_m2',[],'z_m3',[],'z_m4',[]); %vettore dei coefficienti A,B,C,D per le 4 natural frequencies
phi=struct('phi_1',[],'phi_2',[],'phi_3',[],'phi_4',[]); %shape functions corrispondenti alle prime 4 natural frequencies
N=zeros(length(sol_g));

for i=1:4
    
    g=sol_g(i);
    H=[1 0 1 0; 0 1 0 1; -cos(g*L) -sin(g*L) cosh(g*L) sinh(g*L); sin(g*L) -cos(g*L) sinh(g*L) cosh(g*L)];
    N=H(:,1);               %prima colonna della matrice H posta nota poichè scegliamo A=1 arbitrariamente
    H_s=H(:,2:end);         %parte della matrice H che moltiplica per il vettore x_s incognito
    x_s=-(H_s\N);           %-inversa della matrice H restante * N
    campo=sprintf('z_m%d',i);       
    z.(campo)= [1; x_s(1:end)];     %soluzione del vettore z (gli ultimi 3 termini+il primo noto)
    cella=sprintf('phi_%d',i);
    x=sym('x');
    phi.(cella)=matlabFunction(z.(campo)(1)*cos(g*x)+ z.(campo)(2)*sin(g*x)+ z.(campo)(3)*cosh(g*x)+ z.(campo)(4)*sinh(g*x));   %z.(campo)(1) è il primo coeff del vettore z (A B C D) della prima nat freq. e così via)

end

%ho calcolato il mode shape di ciascun modo di vibrare sotto forma di
%funzione nella variabile x

%plot tutti i modi di vibrare in funzione della distanza

x=linspace(0,L,1000);
freq=zeros(length(w_nf),1);     %w_nf è lungo 4, frequenze naturali

figure(3)
for i=1:length(w_nf)
    cella=sprintf('phi_%d',i);
    freq(i)=w_nf(i)/(2*pi);
    y=phi.(cella)(x);
    y_n=y./max(abs(y));     %normalizzazione
    color=['b','r','g','m'];
    subplot(2,2,i)
    plot(x,y_n,color(i),'LineWidth',1.5)
    hold on
    plot(x,zeros(length(x),1),'--k')
    grid on
    xlabel('Distance x from the fixed end [m]')
    ylabel(['Mode shape \phi_{', num2str(i), '}(x)']);
    title(sprintf('mode_{%d} - f_{%d}= %.2f Hz',i,i,freq(i)))
 end


%% Frequency response function

%Calcoliamo le posizione dei nodi per i 4 modi di vibrare, poi utilizzeremo
%queste posizioni per la scelta di xk e xj. Per studiare la reciprocità e
%cosa accade se la forza viene applicata in un nodo.
nodi=struct('n_m1',[],'n_m2',[],'n_m3',[],'n_m4',[]);

%nodes position
for i=1:length(w_nf)
    cella=sprintf('phi_%d',i);
    campo=sprintf('n_m%d',i);
    nodi.(campo)=unique(round(fsolve(phi.(cella),0:0.1:L),4,'decimals'));
    soglia=0.2; %imposto una soglia ed elimino i valori dei nodi sotto di questa
    pos=find((nodi.(campo))<soglia);% trovo le posizioni dei valori sotto soglia e le pongo uguale a zero
    nodi.(campo)(pos)=[];
end

%esempio

output=[0.2 0.6 0.75 0.8]; %[m]       output location
input=[1.2 1.0 0.8 0.3]; %[m]       input location

% output=[0.2 0.6 0.8 0.75]; %[m]       output location
% input=[1.2 1.2 0.8 0.3]; %[m]       input location

% output=[0.2 0.4 0.5 0.7];
% input=[1 0.5 0.1 0.9];
% valutare la G_jk in un nodo e invertire xk e xj per verificare la
% reciprocità, così lo inseriamo nel report.

% modal mass computation
m_i = zeros(length(w_nf),1);

for k=1:length(w_nf)
    cella = sprintf('phi_%d',k);
    m_i(k)= trapz(x,m*(phi.(cella)(x)).^2);
end

%definition of damping ratio
csi = 0.01;
FRF_ex=cell(4,1);
% FRF
for j=1:4
    xk=input(j);
    xj=output(j);

    omega=sym('omega');
    G_jk = 0;
        for i=1:length(w_nf)
            cella = sprintf('phi_%d',i);
            G_jk = matlabFunction (G_jk + ( - (phi.(cella)(xj).*phi.(cella)(xk))/m_i(i) ) / (-(omega)^2 + (2*1j*csi*w_nf(i)*omega) + (w_nf(i))^2)); %mettiamo il meno per la direzione della forza
        end

    FRF_ex{j}=G_jk; % Vettore contenente le quattro G_jk

    freq=linspace(f(1),f(2),5000);
    omega=2*pi*freq;
    figure(j+4)
    subplot(2,1,1)
    semilogy(freq,abs(G_jk(omega)),'b','DisplayName',sprintf('G_{jk}=x_{j}/F_{k}, x_{j}=%.2f m, x_{k}=%.2f m',xj,xk),Linewidth=2)
    xlabel('Frequency [Hz]')
    ylabel('[m/N]')
    grid on
    title('Magnitude')
    legend('show')
    subplot(2,1,2)
    plot(freq,angle(G_jk(omega)),'b',LineWidth=2)
    grid on
    ylabel('[rad]')
    xlabel('Frequency [Hz]')
    title('Phase')
end

%prima G_jk (0.2,1.2)
%seconda G_jk (0.6,1,2)
%Terzo G_jk (0.8,0.8)
%quarta G_jk(0.75,0.3)


%% Model parameters identification - 1ST VIBRATION MODE

n=4;
fmin=3;
fmax=6;
omega_min=2*pi*fmin;
omega_max=2*pi*fmax;
omega_s=linspace(omega_min,omega_max,300)';

G_jk_ex=[FRF_ex{1}(omega_s) FRF_ex{2}(omega_s) FRF_ex{3}(omega_s) FRF_ex{4}(omega_s)]; % matrice dati sperimentali valutando le G_jk ottenute prima nel range di frequenze scelto

syms w_i csi_i
syms A_1  R_L1 R_H1 im_L1 im_H1
syms A_2  R_L2 R_H2 im_L2 im_H2
syms A_3  R_L3 R_H3 im_L3 im_H3
syms A_4  R_L4 R_H4 im_L4 im_H4

G_num_r = cell(4,1);

G_num_r{1}= A_1./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L1+1j*im_L1)./(omega_s).^2)+(R_H1+1j*im_H1);
G_num_r{2}= A_2./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L2+1j*im_L2)./(omega_s).^2)+(R_H2+1j*im_H2);
G_num_r{3}= A_3./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L3+1j*im_L3)./(omega_s).^2)+(R_H3+1j*im_H3);
G_num_r{4}= A_4./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L4+1j*im_L4)./(omega_s).^2)+(R_H4+1j*im_H4);

G_num = [G_num_r{1} G_num_r{2} G_num_r{3} G_num_r{4}];

err_sym = sum(sum(((real(G_jk_ex-G_num)).^2+ (imag(G_jk_ex-G_num)).^2)));

err=matlabFunction(err_sym,'Vars', [w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4]);

% 1) W0 di partenza

w0 = omega_s(151);  %4.5Hz
%ho preso l'elemento 151, supposto che la frequenza w0 sia a metà del range
%scelto per individuare il picco

% 2) Slope method per individuare la csi di partenza

derivata_fase = gradient(angle(FRF_ex{1}(omega_s)),omega_s);
index = find(omega_s==w0);
csi_0 = (-1)/(w0*derivata_fase(index));

% 3) Stima dei MODE SHAPES

A_jk=zeros(n,1);
for k=1:n
    A_jk(k,1)=-imag(FRF_ex{k}(w0))*2*w0^2*csi_0;
end

% Inizio minimizzazione

%Vettore di inizializzazione
%w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4
x0 = [w0,csi_0, A_jk(1), A_jk(2), A_jk(3), A_jk(4), 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0];

% Minimizzazione ai minimi quadrati 

parameters = lsqnonlin(@(x) err(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21),x(22)),x0,[],[],[]); 

omega=linspace(2*pi*4,2*pi*5,200);
G_num_plot=cell(4,1);
freq=linspace(fmin,fmax,300);
figure

for ii=1:n
G_num_plot{ii}= parameters(2+ii)./(-(omega).^2+(1j*2*parameters(2)*parameters(1).*omega)+parameters(1)^2)+((parameters(6+ii)+1j*parameters(10+ii))./(omega).^2)+(parameters(14+ii)+1j*parameters(18+ii));

subplot(2,4,ii)
semilogy(freq,abs(FRF_ex{ii}(omega_s)),'b','DisplayName',sprintf('G^e^x^p_{jk}=x_{j}/F_{k}, x_{j}=%.2f m, x_{k}=%.2f m',output(ii),input(ii)),Linewidth=2)
xlabel('Frequency [Hz]')
ylabel('[m/N]')
ylim([10^(-4), 10^(-1)])
grid on
title('Magnitude')


hold on
semilogy(linspace(4,5,200),abs(G_num_plot{ii}),'or','DisplayName',sprintf('G^n^u^m_%d',ii))
legend('show','Location','south')
subplot(2,4,4+ii)
plot(freq,angle(FRF_ex{ii}(omega_s)),'b',LineWidth=2)
hold on
plot(linspace(4,5,200),angle(G_num_plot{ii}),'or')
grid on
ylabel('[rad]')
xlabel('Frequency [Hz]')
title('Phase')
end
sgtitle('Fitting first vibration mode')


%% Results of identification (1° mode shape)

x=linspace(0,L,1000);
y=phi.phi_1(x);
y_n=y./max(abs(y)); %normalizzazione
Mode_1 = {'Model' ;'Identified'};
NaturalFrequency= {w1/(2*pi);parameters(1)/(2*pi)};
DampingRatio = {csi;csi_0};
% Come mettere il (%) ?
T = table(Mode_1, NaturalFrequency, DampingRatio);

figure
plot(x,y_n,'b','LineWidth',2)
hold on
plot(x,zeros(length(x),1),'--k')
grid on
xlabel('Distance x from the fixed end [m]')
ylabel('Mode shape \phi_{1}(x)')
hold on
point=zeros(n,1);
for ii=1:n
    omega=parameters(1);
    point(ii)=parameters(2+ii)*m_i(1)/phi.phi_1(input(ii)); % restituisce i valori di phi(xj) per le 4 FRF
end

factor = phi.phi_1(output(1))/point(1);

for k=1:n
plot(output(k),point(k)*factor/(max(abs(y))),'o','Color',[0.4 0.7 0.8],'LineWidth',3) %bastava normalizzare questi punti delle phi come quelli analitici
legend('Model','','Identified','location','northeast')
hold on
end
hold on

uitable('Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, 'RowName', T.Properties.RowNames, 'Position', [90 90 240 60],'BackgroundColor',[0.6 0.8 1]);
title('Mode Shape \phi_{1}')

%% Model parameters identification - 2ND VIBRATION MODE

% Stiamo facendo il caso per il secondo picco f2=28.2 Hz

n=4;
fmin=26.2;
fmax=30.2;
omega_min=2*pi*fmin;
omega_max=2*pi*fmax;
omega_s=linspace(omega_min,omega_max,300)';

G_jk_ex=[FRF_ex{1}(omega_s) FRF_ex{2}(omega_s) FRF_ex{3}(omega_s) FRF_ex{4}(omega_s)]; % matrice dati sperimentali valutando le G_jk ottenute prima nel range di frequenze scelto

syms w_i csi_i
syms A_1  R_L1 R_H1 im_L1 im_H1
syms A_2  R_L2 R_H2 im_L2 im_H2
syms A_3  R_L3 R_H3 im_L3 im_H3
syms A_4  R_L4 R_H4 im_L4 im_H4

G_num_r = cell(4,1);

G_num_r{1}= A_1./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L1+1j*im_L1)./(omega_s).^2)+(R_H1+1j*im_H1);
G_num_r{2}= A_2./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L2+1j*im_L2)./(omega_s).^2)+(R_H2+1j*im_H2);
G_num_r{3}= A_3./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L3+1j*im_L3)./(omega_s).^2)+(R_H3+1j*im_H3);
G_num_r{4}= A_4./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L4+1j*im_L4)./(omega_s).^2)+(R_H4+1j*im_H4);

G_num = [G_num_r{1} G_num_r{2} G_num_r{3} G_num_r{4}];

err_sym = sum(sum(((real(G_jk_ex-G_num)).^2+ (imag(G_jk_ex-G_num)).^2)));

err=matlabFunction(err_sym,'Vars', [w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4]);

% 1) W0 di partenza 

w0 = omega_s(151); 
%ho preso l'elemento 151, supposto che la frequenza w0 sia a metà del range
%scelto per individuare il picco

% 2) Slope method per individuare la csi di partenza 

derivata_fase = gradient(angle(FRF_ex{1}(omega_s)),omega_s);
index = find(omega_s==w0);
csi_0 = (-1)/(w0*derivata_fase(index));

% 3) Stima dei MODE SHAPES

A_jk=zeros(n,1);
for k=1:n
    A_jk(k,1)=-imag(FRF_ex{k}(w0))*2*w0^2*csi_0;
end

% Inizio minimizzazione

%Vettore di inizializzazione
%w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4
x0 = [w0,csi_0, A_jk(1), A_jk(2), A_jk(3), A_jk(4), 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0];

% Minimizzazione ai minimi quadrati 

parameters = lsqnonlin(@(x) err(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21),x(22)),x0,[],[],[]); 

omega=linspace(2*pi*(fmin+1),2*pi*(fmax-1),200);
G_num_plot=cell(4,1);
freq=linspace(fmin,fmax,300);
figure

for ii=1:n
G_num_plot{ii}= parameters(2+ii)./(-(omega).^2+(1j*2*parameters(2)*parameters(1).*omega)+parameters(1)^2)+((parameters(6+ii)+1j*parameters(10+ii))./(omega).^2)+(parameters(14+ii)+1j*parameters(18+ii));

subplot(2,4,ii)
semilogy(freq,abs(FRF_ex{ii}(omega_s)),'b','DisplayName',sprintf('G^e^x^p_{jk}=x_{j}/F_{k}, x_{j}=%.2f m, x_{k}=%.2f m',output(ii),input(ii)),Linewidth=2)
xlabel('Frequency [Hz]')
ylabel('[m/N]')
grid on
title('Magnitude')


hold on
semilogy(linspace(fmin+1,fmax-1,200),abs(G_num_plot{ii}),'or','DisplayName',sprintf('G^n^u^m_%d',ii))
legend('show','Location','south')
subplot(2,4,4+ii)
plot(freq,angle(FRF_ex{ii}(omega_s)),'b',LineWidth=2)
hold on
plot(linspace(fmin+1,fmax-1,200),angle(G_num_plot{ii}),'or')
grid on
ylabel('[rad]')
xlabel('Frequency [Hz]')
title('Phase')
end
sgtitle('Fitting second vibration mode')

%% Results of identification (2° mode shape)

x=linspace(0,L,1000);
y=phi.phi_2(x);
y_n=y./max(abs(y)); %normalizzazione
Mode_2 = {'Model' ;'Identified'};
NaturalFrequency= {w2/(2*pi);parameters(1)/(2*pi)};
DampingRatio = {csi;csi_0};
% Come mettere il (%) ?
T = table(Mode_2, NaturalFrequency, DampingRatio);

figure
plot(x,y_n,color(2),'LineWidth',2)
hold on
plot(x,zeros(length(x),1),'--k')
grid on
xlabel('Distance x from the fixed end [m]')
ylabel('Mode shape \phi_{2}(x)')
hold on
point=zeros(n,1);
for ii=1:n
    omega=parameters(1);
    point(ii)=parameters(2+ii)*m_i(2)/phi.phi_2(input(ii)); % restituisce i valori di phi(xj)
end

factor= phi.phi_2(output(1))/point(1);

for k=1:n
plot(output(k),point(k)*factor/max(abs(y)),'o','Color',[0.4 0.7 0.8],'LineWidth',3)
legend('Model','','Identified','location','best')
hold on
end
hold on

uitable('Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, 'RowName', T.Properties.RowNames, 'Position', [90 240 240 60],'BackgroundColor',[0.6 0.8 1]);
title('Mode Shape \phi_{2}')

%% Model parameters identification - 3RD VIBRATION MODE

% Stiamo facendo il caso per il terzo picco f3=79.06 Hz

n=4;
fmin=79.06-2;
fmax=79.06 + 2 ;
omega_min=2*pi*fmin;
omega_max=2*pi*fmax;
omega_s=linspace(omega_min,omega_max,300)';

G_jk_ex=[FRF_ex{1}(omega_s) FRF_ex{2}(omega_s) FRF_ex{3}(omega_s) FRF_ex{4}(omega_s)]; % matrice dati sperimentali valutando le G_jk ottenute prima nel range di frequenze scelto

syms w_i csi_i
syms A_1  R_L1 R_H1 im_L1 im_H1
syms A_2  R_L2 R_H2 im_L2 im_H2
syms A_3  R_L3 R_H3 im_L3 im_H3
syms A_4  R_L4 R_H4 im_L4 im_H4

G_num_r = cell(4,1);

G_num_r{1}= A_1./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L1+1j*im_L1)./(omega_s).^2)+(R_H1+1j*im_H1);
G_num_r{2}= A_2./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L2+1j*im_L2)./(omega_s).^2)+(R_H2+1j*im_H2);
G_num_r{3}= A_3./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L3+1j*im_L3)./(omega_s).^2)+(R_H3+1j*im_H3);
G_num_r{4}= A_4./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L4+1j*im_L4)./(omega_s).^2)+(R_H4+1j*im_H4);

G_num = [G_num_r{1} G_num_r{2} G_num_r{3} G_num_r{4}];

err_sym = sum(sum(((real(G_jk_ex-G_num)).^2+ (imag(G_jk_ex-G_num)).^2)));

err=matlabFunction(err_sym,'Vars', [w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4]);

% 1) W0 di partenza

w0 = omega_s(151);
%ho preso l'elemento 151, supposto che la frequenza w0 sia a metà del range
%scelto per individuare il picco

% 2) Slope method per individuare la csi di partenza 

derivata_fase = gradient(angle(FRF_ex{1}(omega_s)),omega_s);
index = find(omega_s==w0);
csi_0 = (-1)/(w0*derivata_fase(index));

% 3) Stima dei MODE SHAPES

A_jk=zeros(n,1);
for k=1:n
    A_jk(k,1)=-imag(FRF_ex{k}(w0))*2*w0^2*csi_0;
end

% Inizio minimizzazione

%Vettore di inizializzazione
%w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4
x0 = [w0,csi_0, A_jk(1), A_jk(2), A_jk(3), A_jk(4), 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0];

% Minimizzazione ai minimi quadrati 

parameters = lsqnonlin(@(x) err(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21),x(22)),x0,[],[],[]); 

omega=linspace(2*pi*(fmin+1),2*pi*(fmax-1),200);
G_num_plot=cell(4,1);
freq=linspace(fmin,fmax,300);
figure

for ii=1:n
G_num_plot{ii}= parameters(2+ii)./(-(omega).^2+(1j*2*parameters(2)*parameters(1).*omega)+parameters(1)^2)+((parameters(6+ii)+1j*parameters(10+ii))./(omega).^2)+(parameters(14+ii)+1j*parameters(18+ii));

subplot(2,4,ii)
semilogy(freq,abs(FRF_ex{ii}(omega_s)),'b','DisplayName',sprintf('G^e^x^p_{jk}=x_{j}/F_{k}, x_{j}=%.2f m, x_{k}=%.2f m',output(ii),input(ii)),Linewidth=2)
xlabel('Frequency [Hz]')
ylabel('[m/N]')
grid on
title('Magnitude')


hold on
semilogy(linspace(fmin+1,fmax-1,200),abs(G_num_plot{ii}),'or','DisplayName',sprintf('G^n^u^m_%d',ii))
legend('show','Location','south')
subplot(2,4,4+ii)
plot(freq,angle(FRF_ex{ii}(omega_s)),'b',LineWidth=2)
hold on
plot(linspace(fmin+1,fmax-1,200),angle(G_num_plot{ii}),'or')
grid on
ylabel('[rad]')
xlabel('Frequency [Hz]')
title('Phase')
end
sgtitle('Fitting third vibration mode')

%% Results of identification (3° mode shape)

x=linspace(0,L,1000);
y=phi.phi_3(x);
y_n=y./max(abs(y)); %normalizzazione
Mode_3 = {'Model' ;'Identified'};
NaturalFrequency= {w3/(2*pi);parameters(1)/(2*pi)};
DampingRatio = {csi;csi_0};
% Come mettere il (%) ?
T = table(Mode_3, NaturalFrequency, DampingRatio);

figure
plot(x,y_n,color(3),'LineWidth',2)
hold on
plot(x,zeros(length(x),1),'--k')
grid on
xlabel('Distance x from the fixed end [m]')
ylabel('Mode shape \phi_{3}(x)')
hold on
point=zeros(n,1);
for ii=1:n
    omega=parameters(1);
    point(ii)=parameters(2+ii)*m_i(3)/phi.phi_3(input(ii)); % restituisce i valori di phi(xj)
end

factor= phi.phi_3(output(1))/point(1);

for k=1:n
plot(output(k),point(k)*factor/max(abs(y)),'o','Color',[0.4 0.7 0.8],'LineWidth',3)
legend('Model','','Identified','location','northwest')
hold on
end
hold on

uitable('Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, 'RowName', T.Properties.RowNames, 'Position', [230 50 240 60],'BackgroundColor',[0.6 0.8 1]);
title('Mode Shape \phi_{3}')

%% Model parameters identification (4TH VIBRATION MODE)

% Stiamo facendo il caso per il quarto picco f4=154.8 Hz

n=4;
fmin=154.8-2;
fmax=154.8+2;
omega_min=2*pi*fmin;
omega_max=2*pi*fmax;
omega_s=linspace(omega_min,omega_max,300)';

G_jk_ex=[FRF_ex{1}(omega_s) FRF_ex{2}(omega_s) FRF_ex{3}(omega_s) FRF_ex{4}(omega_s)]; % matrice dati sperimentali valutando le G_jk ottenute prima nel range di frequenze scelto

syms w_i csi_i
syms A_1  R_L1 R_H1 im_L1 im_H1
syms A_2  R_L2 R_H2 im_L2 im_H2
syms A_3  R_L3 R_H3 im_L3 im_H3
syms A_4  R_L4 R_H4 im_L4 im_H4

G_num_r = cell(4,1);

G_num_r{1}= A_1./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L1+1j*im_L1)./(omega_s).^2)+(R_H1+1j*im_H1);
G_num_r{2}= A_2./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L2+1j*im_L2)./(omega_s).^2)+(R_H2+1j*im_H2);
G_num_r{3}= A_3./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L3+1j*im_L3)./(omega_s).^2)+(R_H3+1j*im_H3);
G_num_r{4}= A_4./(-(omega_s).^2+(1j*2*csi_i*w_i.*omega_s)+w_i^2)+((R_L4+1j*im_L4)./(omega_s).^2)+(R_H4+1j*im_H4);

G_num = [G_num_r{1} G_num_r{2} G_num_r{3} G_num_r{4}];

err_sym = sum(sum(((real(G_jk_ex-G_num)).^2+ (imag(G_jk_ex-G_num)).^2)));

err=matlabFunction(err_sym,'Vars', [w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4]);

% 1) W0 di partenza

w0 = omega_s(151); 
%ho preso l'elemento 151, supposto che la frequenza w0 sia a metà del range
%scelto per individuare il picco

% 2) Slope method per individuare la csi di partenza 

derivata_fase = gradient(angle(FRF_ex{1}(omega_s)),omega_s);
index = find(omega_s==w0);
csi_0 = (-1)/(w0*derivata_fase(index));

% 3) Stima dei MODE SHAPES

A_jk=zeros(n,1);
for k=1:n
    A_jk(k,1)=-imag(FRF_ex{k}(w0))*2*w0^2*csi_0;
end

% Inizio minimizzazione

%Vettore di inizializzazione
%w_i,csi_i,A_1,A_2,A_3,A_4,R_L1,R_L2,R_L3,R_L4,im_L1,im_L2,im_L3,im_L4,R_H1,R_H2,R_H3,R_H4,im_H1,im_H2,im_H3,im_H4
x0 = [w0,csi_0, A_jk(1), A_jk(2), A_jk(3), A_jk(4), 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0];

% Minimizzazione ai minimi quadrati 

parameters = lsqnonlin(@(x) err(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21),x(22)),x0,[],[],[]); 

omega=linspace(2*pi*(fmin+1),2*pi*(fmax-1),200);
G_num_plot=cell(4,1);
freq=linspace(fmin,fmax,300);
figure

for ii=1:n
G_num_plot{ii}= parameters(2+ii)./(-(omega).^2+(1j*2*parameters(2)*parameters(1).*omega)+parameters(1)^2)+((parameters(6+ii)+1j*parameters(10+ii))./(omega).^2)+(parameters(14+ii)+1j*parameters(18+ii));

subplot(2,4,ii)
semilogy(freq,abs(FRF_ex{ii}(omega_s)),'b','DisplayName',sprintf('G^e^x^p_{jk}=x_{j}/F_{k}, x_{j}=%.2f m, x_{k}=%.2f m',output(ii),input(ii)),Linewidth=2)
xlabel('Frequency [Hz]')
ylabel('[m/N]')
grid on
title('Magnitude')


hold on
semilogy(linspace(fmin+1,fmax-1,200),abs(G_num_plot{ii}),'or','DisplayName',sprintf('G^n^u^m_%d',ii))
legend('show','Location','south')
subplot(2,4,4+ii)
plot(freq,angle(FRF_ex{ii}(omega_s)),'b',LineWidth=2)
hold on
plot(linspace(fmin+1,fmax-1,200),angle(G_num_plot{ii}),'or')
grid on
ylabel('[rad]')
xlabel('Frequency [Hz]')
title('Phase')
end
sgtitle('Fitting fourth vibration mode')

%% Results of identification (4° mode shape)

x=linspace(0,L,1000);
y=phi.phi_4(x);
y_n=y./max(abs(y)); %normalizzazione
Mode_4 = {'Model' ;'Identified'};
NaturalFrequency= {w4/(2*pi);parameters(1)/(2*pi)};
DampingRatio = {csi;csi_0};

T = table(Mode_4, NaturalFrequency, DampingRatio);

figure
plot(x,y_n,color(4),'LineWidth',2)
hold on
plot(x,zeros(length(x),1),'--k')
grid on
xlabel('Distance x from the fixed end [m]')
ylabel('Mode shape \phi_{4}(x)')
hold on
point=zeros(n,1);
for ii=1:n
    omega=parameters(1);
    point(ii)=parameters(2+ii)*m_i(4)/phi.phi_4(input(ii)); % restituisce i valori di phi(xj)
end

factor= phi.phi_4(output(1))/point(1);

for k=1:n
plot(output(k),point(k)*factor/max(abs(y)),'o','Color',[0.4 0.7 0.8],'LineWidth',3)
legend('Model','','Identified','location','northwest')
hold on
end
hold on

uitable('Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, 'RowName', T.Properties.RowNames, 'Position', [300 330 240 60],'BackgroundColor',[0.6 0.8 1]);
title('Mode Shape \phi_{4}')

% FINE