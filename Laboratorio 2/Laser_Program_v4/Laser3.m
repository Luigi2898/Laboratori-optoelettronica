%I programmi implementati per il modello di un laser con 3 popolazioni 
%permettono di fare l'analisi statica e dinamica delle equazioni di bilancio.
%-------------------------------------------------------------------------------------------------------%
function varargout = Laser3(varargin)
%LASER3 M-file for Laser3.fig
%Funzione per la gestione dell'interfaccia di inizializzazione
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                  'gui_Singleton',  gui_Singleton, ...
                  'gui_OpeningFcn', @Laser3_OpeningFcn, ...
                  'gui_OutputFcn',  @Laser3_OutputFcn, ...
                  'gui_LayoutFcn',  [] , ...
                  'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
   [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end
%-------------------------------------------------------------------------------------------------------%
function Laser3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Laser1 (see VARARGIN)
DummyVariable(eventdata);
% Choose default command line output for Laser1
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
%-------------------------------------------------------------------------------------------------------%
function varargout = Laser3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DummyVariable(hObject,eventdata);
% Get default command line output from handles structure
varargout{1} = handles.output;
%-------------------------------------------------------------------------------------------------------%
% DummyVariable to avoid warnings from M-Lint
% It is used to suppress warning from M-Lint for variable not
% used and functions never called. The formers are caused by parameters
% automatically passed by MatLab in callback functions; the latters are
% functions used as callback and therefore passed as a string and
% automatically called by MatLab
function DummyVariable(varargin)
if(0)
   AnalyticAnalyisis_Callback;
   Cavity_A_Callback;
   Material_A_Callback;
   IP_A_Callback;
   RunA_Callback;
   DynamicAnalysis_Callback;
   TypeDynamic_Callback;
   RunD1_Callback;
   IR_Callback;
   RunD2_Callback;
   RunD3_Callback;
   Close_Parameters_Callback;
end
%-------------------------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------------------------%
%Funzione per la chiusura della schermata
function Close_Parameters_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata,handles);

clear all;
close;
%-------------------------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------------------------%
%Analisi statica
function AnalyticAnalyisis_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
%Impostazione parametri interfaccia grafica
set(handles.DynamicAnalysis,'Visible','off'); set(handles.ASD,'Visible','off');
set(handles.PCM_A,'Visible','off'); set(handles.Parameters_A,'Visible','off');
set(handles.Cavity_A,'Visible','off'); set(handles.Material_A,'Visible','off');
set(handles.SP_A,'Visible','off'); set(handles.PSP_A,'Visible','off');
set(handles.IP_A,'Visible','off'); set(handles.T_Ro,'Visible','off');
set(handles.P_Ro,'Visible','off'); set(handles.R_Ro,'Visible','off');
set(handles.P_CER,'Visible','off'); set(handles.T_CER,'Visible','off');
set(handles.R_CER,'Visible','off'); set(handles.DP_D,'Visible','off');
set(handles.Dynamic,'Visible','off'); set(handles.TypeDynamic,'Visible','off');
set(handles.AM,'Visible','off'); set(handles.FM,'Visible','off');
set(handles.RunA,'Visible','off'); set(handles.RunD1,'Visible','off');
set(handles.RunD2,'Visible','off'); set(handles.RunD3,'Visible','off'); 
%Chiamata alla finestra di prompt e impostazione dei suoi parametri
prompt  = {'Wavelenght of Laser3 [nm]'};
dlgTitle   = 'Wavelenght';
lineNo = 1;
def1     = {'1550'};
if isfield(handles, 'answer1')
    def1{1}= num2str(str2num(char(handles.answer1(1))));
end
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
handles.answer1  = inputdlg(prompt,dlgTitle,lineNo,def1,AddOpts);
if isempty(handles.answer1)
   return
end
%Lambda d'operazione  
handles.Lambda=str2double(char(handles.answer1(1)))*1e-7;  % lunghezza d'onda [cm]
%Costanti universali
handles.q=1.6021e-19;  % carica dell' elettrone [C]
handles.c=2.99792e10;   % Velocita della luce nel vuoto [cm/s]
handles.h=6.6261e-34;  % Costante di plack [Js]
guidata(gcbo,handles);
set(handles.ASD,'Visible','on'); set(handles.PCM_A,'Visible','on');
set(handles.Parameters_A,'Visible','on'); set(handles.Cavity_A,'Visible','on');
set(handles.Material_A,'Visible','on');
%-------------------------------------------------------------------------------------------------------%
%Gestione parametri cavità
function Cavity_A_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
%Funzioni di impostazione parametri interfaccia grafica
handles.OPR=[];
set(handles.IP_A,'Visible','off'); set(handles.PSP_A,'Visible','off');
set(handles.SP_A,'Visible','off'); set(handles.P_Ro,'Visible','off'); 
set(handles.T_Ro,'Visible','off'); set(handles.R_Ro,'Visible','off');
set(handles.RunA,'Visible','off');
%Chiamata alla finestra di prompt e impostazione dei suoi parametri
prompt  = {'Depth [nm]','Width [um]','Active Region Length [um]',...
           'Passive Region Length [um]',...
           'Separate Confinement Thickness [um]',...
           '\Gamma _{XY} --> Value of Gama total in the well ',...
           '\Gammab _{XY} --> Value of Gama total in the barrier ',...
           '\beta _{sp} --> Spontaneous Emission Factor',...
           'R1 --> Reflectivity Minimum',...
           'R2 --> Reflectivity Maximum','\eta_g --> Group Index'};
dlgTitle='Cavity Parameters';
lineNo=1;
def2={'10','2','250','0','0.3','0.032','0','[1e-4]','[0.32]','[0.32]','4.2'};
if isfield(handles, 'answer2')
   def2=handles.answer2;
end
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
handles.answer2  = inputdlg(prompt,dlgTitle,lineNo,def2,AddOpts);
if isempty(handles.answer2)
   return
end
%Costanti cavità
handles.d=str2double(char(handles.answer2(1)))*1e-7;   % Spessore complesivo del well o dell'attivo [cm]
handles.w=str2double(char(handles.answer2(2)))*1e-4;   % Larghezza della cavita' [cm]
handles.La=str2double(char(handles.answer2(3)))*1e-4;  % Lunghezza attiva della cavita' [cm]
handles.Lp=str2double(char(handles.answer2(4)))*1e-4;  % Lunghezza passiva della cavita' [cm]
handles.Lsch=str2double(char(handles.answer2(5)))*1e-4;
handles.Gxy=str2double(char(handles.answer2(6))); % Valori di gamma complessivo di tutti well
handles.Gxyb=str2double(char(handles.answer2(7)));   % Valori del gamma complessivo nella barriera 
handles.B1=str2num(char(handles.answer2(8)));     % Fattore di emissione spontanea
handles.B1_Originale=handles.B1;
handles.R1=str2num(char(handles.answer2(9)));     % Riflettivita minima
handles.R1_Originale=handles.R1;
handles.R2=str2num(char(handles.answer2(10)));    % Riflettivita masssima
handles.R2_Originale=handles.R2;
handles.ng=str2double(char(handles.answer2(11))); % indice di gruppo
%Calcolo potenza
for m=1:length(handles.R1)
    for n=1:length(handles.R2)
        handles.F1(m,n) = (1-handles.R1(m)^2)/((1-handles.R1(m)^2)+...   
                           handles.R1(m)/handles.R2(n)*(1-handles.R2(n)^2));   
        handles.OPR(m,n) = handles.F1(m,n)./(1-handles.F1(m,n));  %Output Facest Power Ratio
    end
end
handles.OPR = handles.OPR(:);
set(handles.P_Ro,'Visible','on'); set(handles.T_Ro,'Visible','on');
set(handles.R_Ro,'Visible','on'); set(handles.R_Ro,'String',num2str(handles.OPR));
%Memorizza i dati nell'interfaccia
guidata(gcbo,handles);
if(prod(size(handles.OPR))<1)
    %Display report in second figure
    fig2 = figure('menubar','none','name','Output Facets Power Ratio');
    specs=flipud(handles.R(:));
    m=length(specs);    % Defining layout constants - change to adjust 'look and feel'The names of the tests 
    TestNames = handles.R; % Number of test columns 
    NumTests = length(TestNames);
    NumRows = m+1;      % Total number of rows - header (1) + number of results (m)
    TopMargin = 0.05;   % Margin between top of figure and title row 
    BotMargin = 0.20;   % Margin between last test row and bottom of figure 
    LftMargin = 0.03;   % Margin between left side of figure and Computer Name 
    RgtMargin = 0.03;   % Margin between last test column and right side of figure 
    CNWidth = 0.10;     % Width of Computer Name column 
    MidMargin = 0.03;   % Margin between Computer Name column and first test column
    HBetween = 0.025;   % Distance between two rows of tests
    WBetween = 0.015;   % Distance between two columns of tests
    % Width of each test column
    TestWidth = (1-LftMargin-CNWidth-MidMargin-RgtMargin-(NumTests-1)*WBetween)/NumTests;
    % Height of each test row
    RowHeight = (1-TopMargin-(NumRows-1)*HBetween-BotMargin)/NumRows;
    % Beginning of first test column
    BeginTestCol = LftMargin+CNWidth+MidMargin;
    % Retrieve the background color for the figure
    bc = get(fig2,'Color');
    % Create headers
    % Computer Name column header
    uicontrol('Style', 'text', 'Units', 'normalized', ...
        'Position', [LftMargin 1-TopMargin-RowHeight CNWidth RowHeight],...
        'String', 'R1\R2', 'BackgroundColor', bc, 'FontWeight','bold');
    % Test name column header
    for k=1:NumTests
        uicontrol('Style', 'text', 'Units', 'normalized', 'Position', ...
            [BeginTestCol+(k-1)*(WBetween+TestWidth) 1-TopMargin-RowHeight TestWidth RowHeight],...
            'String',  num2str(TestNames(k)), 'BackgroundColor', bc, ...
            'FontWeight', 'bold');
    end
    % For each computer
    for k=1:NumRows-1
        VertPos = 1-TopMargin-k*(RowHeight+HBetween)-RowHeight;
        % Computer Name row header
        uicontrol('Style', 'text', 'Units', 'normalized', 'Position', ...
            [LftMargin VertPos CNWidth RowHeight],...
            'String', num2str(specs(NumRows-k)), 'BackgroundColor', bc,...
            'HorizontalAlignment', 'left','FontWeight', 'bold');
        % Test results for that computer
        for n=1:NumTests
            uicontrol('Style', 'text', 'Units', 'normalized', 'Position', ...
                [BeginTestCol+(n-1)*(WBetween+TestWidth) VertPos TestWidth RowHeight],...
                'String', sprintf('%.4f',handles.OPR(NumRows-k, n)), 'BackgroundColor', bc);
        end
    end
end    
%-------------------------------------------------------------------------------------------------------%
%Gestione parametri materiali
%Impostazione dei parametri relativi ai materiali per cui si vuole effettuare l'analisi
function Material_A_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
%Funzioni di impostazione parametri interfaccia grafica
set(handles.RunA,'Visible','off');
%Chiamata alla finestra di prompt e impostazione dei suoi parametri
prompt  = {'Differential Gain [cm^2] --> a0',...
           '\eta_i --> Quantum Efficiency',...
           'N_{wo} [cm^{-3}] --> Transparency Carrier Density',...
           '\epsilon_1 [cm^3] --> Gain Compresssion Factor',...
           '\epsilon_2 [cm^6] --> Gain Compresssion Factor',...
           'A [s^{-1}]  -> Monomolecular Recombination Coefficient',...
           'B [cm^3/s] --> Bimolecular Recombination Coefficient',...
           'C [cm^6/s] --> Auger Recombination Coefficient'};
dlgTitle   = 'Materials Parameters';
lineNo = 1;
def3     = {'[5.34e-16]','0.8','1.8e18','1.5e-16',...
            '0e-32','3.57e8','0.8e-10','3.5e-30'};
if isfield(handles, 'answer3')
   def3=handles.answer3;
end
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
handles.answer3  = inputdlg(prompt,dlgTitle,lineNo,def3,AddOpts);
if isempty(handles.answer3)
   return
end
prompt  = {'\tau_d [ns] --> SCH diffusin time constant [ns]',...
           '\tau_g [ns] --> Gateway state lifetime [ns]',...
           '\tau_{ns} [ns] --> Non-radiative carrier lifetime SCH [ns]',...
           '\tau_{ng} [ns] --> Non-radiative carrier lifetime gateway[ns]',...
           '\tau_{esc} [ns] --> Ambipolar escape time [ns]',...
           '\tau_{cap} [ns] --> Ambipolar capture time [ns]',...
           '\alpha_{a} [cm^{-1}] --> Internal Losses in the Active region',...
           '\alpha_{p} [cm^{-1}] --> Internal Losses in the Passive region',...
           '\alpha_{H} LEF --> Linewidth Enhancement Factor'};
dlgTitle   = 'Materials Parameters';
lineNo = 1;
def4     = {'45e-3','1e-3','1','1','5.6e-3','1e-3','[5]','[0]','[4]'};
if isfield(handles, 'answer4')
   def4=handles.answer4;
end
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
handles.answer4  = inputdlg(prompt,dlgTitle,lineNo,def4,AddOpts);
if isempty(handles.answer4)
   return
end
%Costanti del materiale  InGaAs/GaAs
handles.a0=str2num(char(handles.answer3(1)));     % Guadagno differenziale [cm^2]
handles.a0_Originale=handles.a0; 
handles.ni=str2double(char(handles.answer3(2)));  % Efficienza quantica
handles.Nwo=str2double(char(handles.answer3(3))); % Transparency carrier density [cm^-3]
handles.Ep=str2num(char(handles.answer3(4)));  % Gain compression coefficient for photon density [cm^3]
handles.Ep_Originale=handles.Ep;
handles.Ep2=str2double(char(handles.answer3(5))); % 
handles.A=str2double(char(handles.answer3(6)));  % [s^-1]
handles.B=str2double(char(handles.answer3(7)));  % Recombinazione monomoleculare [cm^3/s]
handles.C=str2double(char(handles.answer3(8)));  % Recombinazione bimoleculare [cm^6/s]
handles.Taud=str2double(char(handles.answer4(1)))*1e-9;% SCH diffusin time constant [s]
handles.Taug=str2double(char(handles.answer4(2)))*1e-9;% Gateway state lifetime [s]
handles.Tauns=str2double(char(handles.answer4(3)))*1e-9;% Non-radiative carrier lifetime SCH [s]
handles.Taung=str2double(char(handles.answer4(4)))*1e-9;% Non-radiative carrier lifetime gateway[s]
handles.Tauesc=str2double(char(handles.answer4(5)))*1e-9;% Ambipolar escape time [s]
handles.Taucap=str2double(char(handles.answer4(6)))*1e-9;% Electron capture time [s]
handles.Alf_a=str2num(char(handles.answer4(7)));% Perdite nella parte attiva [cm^-1]
handles.Alf_a_Originale=handles.Alf_a;
handles.Alf_p=str2num(char(handles.answer4(8)));% Perdite nella parte passiva [cm^-1]
handles.Alf_p_Originale=handles.Alf_p;
handles.Alf_H=str2num(char(handles.answer4(9)));% LEF, Efetto allargamento di riga
handles.Alf_H_Originale=handles.Alf_H;
%Memorizza i dati nell'interfaccia
guidata(gcbo,handles);
%Rende visibili i tasti utilizzabili da adesso
set(handles.PSP_A,'Visible','on'); set(handles.SP_A,'Visible','on');
set(handles.IP_A,'Visible','on');
%-------------------------------------------------------------------------------------------------------%
%Gestione Potenza di Uscita
%Si specifica il range della Pout e i punti di analisi
function IP_A_Callback(hObject, eventdata, handles)

%Chiamata alla finestra di prompt e impostazione dei suoi parametri
prompt  = {'Power min [mw]','Power max [mw]','Number of Points'};
dlgTitle   = 'Output Power Range';
lineNo = 1;
def5    = {'0','20','100'};
if isfield(handles, 'answer5')
   def5=handles.answer5;
end
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
handles.answer5  = inputdlg(prompt,dlgTitle,lineNo,def5,AddOpts);
if isempty(handles.answer5)
   return
end
%Parametri di potenza
handles.Poutmin = str2double(char(handles.answer5(1)))*1e-3;    %Potenza minima
handles.Poutmax = str2double(char(handles.answer5(2)))*1e-3;    %Potenza Massima
handles.NumPoints = fix(str2double(char(handles.answer5(3))));  %Numero di punti di simulazione
%Limitazione del valore di potenza minimo e del numero di punti della simulazione
if handles.Poutmin < 1e-12
   handles.Poutmin = 1e-12;
end
if handles.NumPoints < 5
   handles.NumPoints = 10;
end
handles.Pout=linspace(handles.Poutmin,handles.Poutmax,handles.NumPoints);
handles.dnb=0;
handles.dns=0;
handles.dng=0;
%Calcolo delle variabili necessarie per le equazioni finali
handles.Vg=handles.c/handles.ng;                % Velocita di Gruppo [m/s]
handles.v=handles.c/handles.Lambda;             % Frequenza del laser1 [Hz]
handles.Vw=handles.d*handles.w*handles.La;      % Volume della cavita' [cm^3]
handles.Vsch=handles.w*handles.La*handles.Lsch;  % SCH Volume [cm^3]
handles.L=handles.La+handles.Lp;                % Lunghezza della cavita'
handles.Gamw=handles.Gxy*(handles.La/handles.L); % Gama, Fattore di confinamento elettrone - fotoni
handles.Gamb=handles.Gxyb*(handles.La/handles.L);% Confinamento nella barriera
handles.Gams=0; % Confinement factor in total SCH
handles.Gamg=handles.Gamw;% Confinement factor in total gateway
%Calcolo delle equazioni relative alle perdite interne, del vetro e altre
handles.F = (1-handles.R1)/((1-handles.R1)+...
            (1-handles.R2)*sqrt(handles.R1/handles.R2));    % Rapporto Pout/Pouttotale
handles.R = sqrt(handles.R1*handles.R2);                    % Reflettivita' sqrt(R1R2)
handles.Alf_m = (1/handles.L)*log(1/handles.R);             % Perdite dello specchio [m^-1]
handles.Alf_i = (handles.Alf_a*handles.La+handles.Alf_p*handles.Lp)/...
                 handles.L;                                 % Perdite medie interne della cavita [m^-1]
handles.n0 = handles.F*(handles.Alf_m/...
            (handles.Alf_m+handles.Alf_i));                  % Eficenzia Ottica
handles.gth = (handles.Alf_m+handles.Alf_i)/handles.Gamw;    % Guadragno del materiale alla soglia [m^-1]
handles.Taup = 1/(handles.Vg*(handles.Alf_i+handles.Alf_m));% Tempo di vita dei fotoni [s]
handles.Nwnth = handles.Nwo+handles.gth/handles.a0;
handles.Nwth = handles.Nwnth*handles.Vw;
handles.Taunw = 1/(handles.A+handles.B*handles.Nwth+handles.C*handles.Nwth^2);
handles.Ngth = handles.Taucap*handles.Nwth*(1/handles.Tauesc+1/handles.Taunw);
handles.Nschth = handles.Taud*(handles.Ngth.*(1/handles.Taucap+...
                 1/handles.Taug+1/handles.Taung)-handles.Nwth./handles.Tauesc);
handles.Ith = handles.q/handles.ni*(handles.Nschth*(1/handles.Tauns+1/handles.Taud)-...
              handles.Ngth/handles.Taug); %Current of Threshold
guidata(gcbo,handles);
set(handles.RunA,'Visible','on');
%-------------------------------------------------------------------------------------------------------%
%Funzione per la generazione dei grafici relativi l'analisi statica
%I grafici ottenuti riguardano P-I, I-Carrier Density, I-Linewidth
function RunA_Callback(hObject, eventdata, handles)

DummyVariable(eventdata);
handles.val=0; handles.case=0;
set(get(hObject,'parent'),'Pointer','watch');
handles.Nwnop=[];handles.Nwop=[];handles.Taunw=[];handles.Ngop=[];
handles.Nschop=[];handles.I=[];handles.Av=[];
vettorescansione=''; handles.quantitascansione='';
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos=[5,35,scnsize(3)-8,scnsize(4)-115];
handles.stan=figure('Name','Static analysis','NumberTitle','off','Position',pos);
handles.S = handles.Pout*handles.Taup./(handles.n0*handles.h*handles.v); % Total photon number 
handles.go = handles.gth.*(1+handles.Ep*handles.Gamw*...
             handles.S+handles.Ep2*(handles.Gamw*handles.S).^2);    % Saturated Gain
if(length(handles.B1)>=1)
    vettorescansione=handles.B1;handles.quantitascansione='\beta';
elseif(length(handles.R1)>1)
    vettorescansione=handles.R1;handles.quantitascansione='R_1';
elseif(length(handles.R2)>1)
    vettorescansione=handles.R2;handles.quantitascansione='R_2';
elseif(length(handles.a0)>1)
    vettorescansione=handles.a0;handles.quantitascansione='a0';
elseif(length(handles.Alf_a)>1)
    vettorescansione=handles.Alf_a;handles.quantitascansione='\alpha_a';
elseif(length(handles.Alf_p)>1)
    vettorescansione=handles.Alf_p;handles.quantitascansione='\alpha_p';
elseif(length(handles.Alf_H)>1)
    vettorescansione=handles.Alf_H;handles.quantitascansione='\alpha_H';    
end
if strcmp(vettorescansione,'')==0 && strcmp(handles.quantitascansione,'')==0 
    BB1 = length(vettorescansione);
    cicli = length(handles.Pout); % Numeri di cicli
    for j=1:1:BB1
        handles.I=[];handles.Nsch=[];handles.Ng=[];handles.Nw=[];handles.Poutn=[];
        switch(handles.quantitascansione)
            case {'\beta'}
                handles.B1=vettorescansione(j);
                handles.ValorePerlegenda=handles.B1;
            case {'R_1'}
                handles.R1=vettorescansione(j);
                handles.ValorePerlegenda=handles.R1;
            case {'R_2'}
                handles.R2=vettorescansione(j);
                handles.ValorePerlegenda=handles.R2;        
            case {'a0'}
                handles.a0=vettorescansione(j);
                handles.ValorePerlegenda=handles.a0;
            case {'\alpha_a'}
                handles.Alf_a=vettorescansione(j);
                handles.ValorePerlegenda=handles.Alf_a;          
            case {'\alpha_p'}
                handles.Alf_p=vettorescansione(j);
                handles.ValorePerlegenda=handles.Alf_p;
            case {'\alpha_H'}
                handles.Alf_H=vettorescansione(j);
                handles.ValorePerlegenda=handles.Alf_H;          
        end
        %Calcolo dei valori necessari per la creazione dei grafici nei vari punti
        for i=1:1:cicli
            handles.Nwnop(i) = handles.Nwo+handles.go(i)./handles.a0; % Nwop denied
            handles.Nwop(i) = handles.Nwnop(i)*handles.Vw;            % Nwop
            handles.Taunw(i) = 1/(handles.A+handles.B*handles.Nwop(i)+handles.C*handles.Nwop(i)^2);
            handles.Ngop(i) = handles.Taucap*(handles.Nwop(i).*(1/handles.Tauesc+1/...
                              handles.Taunw(i))+(handles.Vg*handles.go(i)*handles.Gamw*...
                              (handles.S(i)/(1+handles.Ep*handles.Gamw*handles.S(i)+...
                              handles.Ep2*(handles.Gamw*handles.S(i)).^2))));
            handles.Nschop(i) = handles.Taud*(handles.Ngop(i).*(1/handles.Taucap+1/...
                                handles.Taug+1/handles.Taung)-handles.Nwop(i)./handles.Tauesc);
            handles.I(i) = handles.q/handles.ni*(handles.Nschop(i).*(1/handles.Tauns+1/...
                           handles.Taud)-handles.Ngop(i)./handles.Taug);
            handles.Av(i) = (handles.Alf_H.*handles.Vg.*handles.Gamw*handles.a0*...
                            handles.Nwop(i)/(4*handles.Vw*pi))+(handles.v/handles.ng)*...
                            (handles.Gams*handles.dns*handles.Nschop(i)/handles.Vsch+...
                            handles.Gamg*handles.dng*handles.Ngop(i)/handles.Vw); 
        end
        %creazione vettori con valore iniziale 0 e finale 50
        handles.Nw = linspace(0,handles.Nwth,50);
        handles.Ng = linspace(0,handles.Ngth,50);
        handles.Nsch = linspace(0,handles.Nschth,50);
        handles.Is = linspace(0,handles.Ith,50);
        handles.I = [handles.Is handles.I];
        handles.Poutn = [zeros(1,50) handles.Pout];
        handles.Nw = [handles.Nw handles.Nwop];
        handles.Ng = [handles.Ng handles.Ngop];
        handles.Nsch = [handles.Nsch handles.Nschop];
        %Creazione dei tre grafici con relative legende
        if(nargout==0)
            GraficiA(j,handles);
        end
    end
end
%ripristino i valori del vettore
switch(handles.quantitascansione)
    case {'\beta'}
        handles.B1=vettorescansione;
    case {'R_1'}
        handles.R1=vettorescansione;
    case {'R_2'}
        handles.R2=vettorescansione;
    case {'a0'}
        handles.a0=vettorescansione;
    case {'\alpha_a'}
        handles.Alf_a=vettorescansione;
    case {'\alpha_p'}
        handles.Alf_p=vettorescansione;
    case {'\alpha_H'}
        handles.Alf_H=vettorescansione;
end
guidata(gcbo,handles);
set(handles.DynamicAnalysis,'Visible','on');
set(get(hObject,'parent'),'Pointer','arrow');
%-------------------------------------------------------------------------------------------------------%
%Funzione che definisce i colori per le diverse curve dei grafici
%Chiamata dopo la creazione di ogni grafico
function gradeplots(h,groupingfactor)

if nargin<1
   h=gca;
end
if nargin<2
   groupingfactor=1;
end
lines=findobj(h,'type','line');
c=colormap;
indici=fix(linspace(1,size(c,1),length(lines)/groupingfactor));
c=c(indici,:);
for m=1:length(lines)
   set(lines(m),'color',c(ceil(m/groupingfactor),:));
end
%-------------------------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------------------------%
%Analisi Dinamica
function DynamicAnalysis_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
%Impostazione parametri interfaccia grafica
set(handles.ASD,'Visible','off'); set(handles.PCM_A,'Visible','off');
set(handles.Parameters_A,'Visible','off'); set(handles.Cavity_A,'Visible','off');
set(handles.Material_A,'Visible','off'); set(handles.SP_A,'Visible','off');
set(handles.PSP_A,'Visible','off');  set(handles.IP_A,'Visible','off');  
set(handles.DP_D,'Visible','off'); set(handles.Dynamic,'Visible','off');
set(handles.TypeDynamic,'Visible','off'); set(handles.AM,'Visible','off');
set(handles.FM,'Visible','off'); set(handles.R_CER,'Visible','off');
set(handles.T_CER,'Visible','off'); set(handles.P_CER,'Visible','off');
set(handles.RunA,'Visible','off'); set(handles.RunD1,'Visible','off');
set(handles.RunD2,'Visible','off'); set(handles.RunD3,'Visible','off'); 
%Associazione valori originali utilizzati durante l'analisi statica
StringaRichiesta='';
defaultvalue='none';
if(length(handles.B1_Originale)>1)
    StringaRichiesta='\beta _{sp} --> Spontaneous Emission';
    defaultvalue=handles.B1_Originale(1);
elseif(length(handles.R1_Originale)>1)
    StringaRichiesta='R1 --> Reflectivity Minimum';
    defaultvalue=handles.R1_Originale(1);
elseif(length(handles.R2_Originale)>1)
    StringaRichiesta='R2 --> Reflectivity Maximum';
    defaultvalue=handles.R2_Originale(1);
elseif(length(handles.a0_Originale(1))>1)
    StringaRichiesta='Differential Gain [cm^2] --> a0(N-N_{tr})';
    defaultvalue=handles.a0_Originale(1);
elseif(length(handles.Alf_a_Originale)>1)
    StringaRichiesta='\alpha_a [cm^{-1}] --> Internal Losses in the Active Region';
    defaultvalue=handles.Alf_a_Originale(1);
elseif(length(handles.Alf_p_Originale)>1)
    StringaRichiesta='\alpha_p [cm^{-1}] --> Internal Losses in the Passive Region';
    defaultvalue=handles.Alf_p_Originale(1);
elseif(length(handles.Alf_H_Originale)>1)
    StringaRichiesta='\alpha_H --> Linewidth Enhancement Factor';
    defaultvalue=handles.Alf_H_Originale(1);
end
%Chiamata alla finestra di prompt e impostazione dei suoi parametri
prompt  = {'Step Start Time [ns]','End Time [ns]'};
dlgTitle = 'Simulation Parameters ';
lineNo = 1;
def6 = {'2','15'};
if strcmp(defaultvalue,'none')==0
    prompt{end+1} = StringaRichiesta;
    def6{end+1} = num2str(defaultvalue);
end
if isfield(handles, 'answer6')
   def6=handles.answer6;
end
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
handles.answer6  = inputdlg(prompt,dlgTitle,lineNo,def6,AddOpts);
if isempty(handles.answer6)
   return
end
%Impostazione tempi di analisi
handles.Tstart = 0;                                    %Tempo di partenza
handles.Tst = str2double(char(handles.answer6(1)))*1e-9;    %Valore dello step  
handles.Tend = str2double(char(handles.answer6(2)))*1e-9;   %Tempo di fine
if handles.Tst <= 1e-9
   handles.Tst = 1e-9;
end
%Calcolo delle variabili necessarie per le equazioni finali
handles.Vg=handles.c/handles.ng;                % Velocita di Gruppo [m/s]
handles.v=handles.c/handles.Lambda;             % Frequenza del laser1 [Hz]
handles.Vw=handles.d*handles.w*handles.La;      % Volume della cavita' [cm^3]
handles.Vsch=handles.w*handles.La*handles.Lsch;  % SCH Volume [cm^3]
handles.L=handles.La+handles.Lp;                % Lunghezza della cavita'
handles.Gamw=handles.Gxy*(handles.La/handles.L); % Gama, Fattore di confinamento elettrone - fotoni
handles.Gamb=handles.Gxyb*(handles.La/handles.L);% Confinamento nella barriera
handles.Gams=0; % Confinement factor in total SCH
handles.Gamg=handles.Gamw;% Confinement factor in total gateway
%Calcolo delle equazioni relative alle perdite interne, del vetro e altre
handles.Alf_m = (1/handles.L)*log(1/handles.R);             % Perdite dello specchio [m^-1]
handles.Alf_i = (handles.Alf_a*handles.La+handles.Alf_p*handles.Lp)/...
                 handles.L;                                 % Perdite medie interne della cavita [m^-1]
handles.n0 = handles.F*(handles.Alf_m/...
            (handles.Alf_m+handles.Alf_i));                  % Eficenzia Ottica
handles.gth = (handles.Alf_m+handles.Alf_i)/handles.Gamw;    % Guadragno del materiale alla soglia [m^-1]
handles.Taup = 1/(handles.Vg*(handles.Alf_i+handles.Alf_m));% Tempo di vita dei fotoni [s]
handles.Nwnth = handles.Nwo+handles.gth/handles.a0;
handles.Nwth = handles.Nwnth*handles.Vw;
handles.Taunw = 1/(handles.A+handles.B*handles.Nwth+handles.C*handles.Nwth^2);
handles.Ngth = handles.Taucap*handles.Nwth*(1/handles.Tauesc+1/handles.Taunw);
handles.Nschth = handles.Taud*(handles.Ngth.*(1/handles.Taucap+...
                 1/handles.Taug+1/handles.Taung)-handles.Nwth./handles.Tauesc);
handles.Ith=[];
[Ith]=GetSoglia(handles.I,handles.Poutn);
handles.Ith=Ith;
handles.Rsp=handles.Nwth/handles.Taunw;               % Spontaneos emission
handles.Espont=handles.B1*handles.Gamw*handles.Rsp;
if strcmp(defaultvalue,'none')==0
    switch (StringaRichiesta)
        case '\beta _{sp} --> Spontaneous Emission'
            handles.B1 = str2double(char(handles.answer5(3)));
        case 'R1 --> Reflectivity Minimum'
            handles.R1 = str2double(char(handles.answer5(3)));
        case 'R2 --> Reflectivity Maximum'
            handles.R2 = str2double(char(handles.answer5(3)));
        case 'Differential Gain [cm^2] --> a0(N-N_{tr})'
    	    handles.a0 = str2double(char(handles.answer5(3)));
        case '\alpha_a [cm^{-1}] --> Internal Losses in the Active Region'
    	    handles.Alf_a = str2double(char(handles.answer5(3)));
        case '\alpha_p [cm^{-1}] --> Internal Losses in the Passive Region'
    	    handles.Alf_p = str2double(char(handles.answer5(3)));
        case '\alpha_H --> Linewidth Enhancement Factor'
    	    handles.Alf_H = str2double(char(handles.answer5(3)));   
    end
end
guidata(gcbo,handles);
set(handles.ASD,'Visible','on'); set(handles.DP_D,'Visible','on');
set(handles.Dynamic,'Visible','on'); set(handles.TypeDynamic,'Visible','on');
set(handles.DP_D,'Position',get(handles.PCM_A,'Position'));
set(handles.Dynamic,'Position',get(handles.Parameters_A,'Position'));
set(handles.TypeDynamic,'Position',get(handles.Cavity_A,'Position'));
%-------------------------------------------------------------------------------------------------------%
%Scelta del tipo di analisi dinamica
function TypeDynamic_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
switch get(handles.TypeDynamic,'Value')
    case 1
        set(handles.AM,'Visible','off'); set(handles.FM,'Visible','off');
        set(handles.RunA,'Visible','off'); set(handles.RunD1,'Visible','off');
        set(handles.RunD2,'Visible','off'); set(handles.RunD3,'Visible','off');
        set(handles.T_CER,'Visible','off'); set(handles.P_CER,'Visible','off');
        set(handles.R_CER,'Visible','off');
        handles.case=1;
    case 2    %Analisi Electro-Optic
       %Impostazione parametri interfaccia grafica
       set(handles.AM,'Visible','off'); set(handles.FM,'Visible','off');
       set(handles.RunA,'Visible','off'); set(handles.RunD2,'Visible','off');
       set(handles.RunD3,'Visible','off'); set(handles.T_CER,'Visible','off');
       set(handles.P_CER,'Visible','off'); set(handles.R_CER,'Visible','off');
       %Chiamata alla finestra di prompt e impostazione dei suoi parametri
       prompt = {'Bias Current  / Ith'};
       dlgTitle = 'Electro-Optic Parameter';
       lineNo = 1;
       def7 = {'2'};
       if isfield(handles, 'answer7')
           def7=handles.answer7;
       end
       AddOpts.Resize = 'on';
       AddOpts.WindowStyle = 'normal';
       AddOpts.Interpreter = 'tex';
       handles.answer7 = inputdlg(prompt,dlgTitle,lineNo,def7,AddOpts);
       if isempty(handles.answer7)
           return
       end
       handles.Bias = str2double(char(handles.answer7(1)));
       handles.Imincs = handles.Bias*handles.Ith;
       handles.Imax = handles.Imincs*0.01;                                 
       handles.Imaxcs=handles.Imincs+handles.Imax;
       handles.title=1;
       guidata(gcbo,handles);
       set(handles.RunD1,'Visible','on'); 
       set(handles.RunD1,'Position',get(handles.RunA,'Position'));
   case 3   %Analisi Switching On
       %Impostazione parametri interfaccia grafica
       set(handles.AM,'Visible','off'); set(handles.FM,'Visible','off');
       set(handles.RunA,'Visible','off'); set(handles.RunD1,'Visible','off');
       set(handles.RunD3,'Visible','off'); set(handles.T_CER,'Visible','off');
       set(handles.P_CER,'Visible','off'); set(handles.R_CER,'Visible','off');       
       %Chiamata alla finestra di prompt e impostazione dei suoi parametri
       prompt = {'Minimun Current / Ith ','Maximun Current / Ith'};
       dlgTitle = 'Switching On Parameter';
       lineNo = 1;
       def8 = {'0.2','2'};
       if isfield(handles, 'answer8')
           def8=handles.answer8;
       end
       AddOpts.Resize = 'on';
       AddOpts.WindowStyle = 'normal';
       AddOpts.Interpreter = 'tex';
       handles.answer8 = inputdlg(prompt,dlgTitle,lineNo,def8,AddOpts);
       if isempty(handles.answer8)
           return
       end
       handles.Imin = str2double(char(handles.answer8(1)));
       handles.Imax = str2double(char(handles.answer8(2)));
       handles.Imincs=handles.Imin*handles.Ith;
       handles.Imaxcs=handles.Imax*handles.Ith;
       handles.title=2;
       guidata(gcbo,handles);
       set(handles.RunD2,'Visible','on'); set(handles.RunD2,'Position',get(handles.RunA,'Position'));
   case 4   %Analisi PRBS Signal
       %Impostazione parametri interfaccia grafica
       set(handles.RunD1,'Visible','off'); set(handles.RunD2,'Visible','off');
       set(handles.AM,'Visible','off'); set(handles.FM,'Visible','off');
       %Chiamata alla finestra di prompt e impostazione dei suoi parametri
       prompt = {'Minimun Current / Ith','Maximun Current / Ith',...
               'Bit Rate [Gbit/s]','Bit Number --> (2"n -1)'};
       dlgTitle = 'PRBS Signal Parameter';
       lineNo = 1;
       def9 = {'1.5','2.5','1','5'};     
       if isfield(handles, 'answer9')
           def9=handles.answer9;
       end
       AddOpts.Resize = 'on';
       AddOpts.WindowStyle = 'normal';
       AddOpts.Interpreter = 'tex';
       handles.answer9 = inputdlg(prompt,dlgTitle,lineNo,def9,AddOpts);
       if isempty(handles.answer9)
           return
       end
       handles.Imin = str2double(char(handles.answer9(1)));
       handles.Imax = str2double(char(handles.answer9(2)));
       handles.Rb = str2double(char(handles.answer9(3)))*1e9;
       handles.Nbit = str2double(char(handles.answer9(4)));
       handles.Imincs=handles.Imin*handles.Ith;
       handles.Imaxcs=handles.Imax*handles.Ith;
       handles.title=3;
       guidata(gcbo,handles);
       set(handles.RunD3,'Visible','on'); 
       set(handles.RunD3,'Position',get(handles.RunA,'Position'));
end
if handles.val==1
    delete(handles.s11);delete(handles.s12);delete(handles.s21);
    delete(handles.s22);delete(handles.s31);delete(handles.s32);
    delete(handles.s41);delete(handles.s42);delete(handles.t11);
    delete(handles.t12);delete(handles.t21);delete(handles.t22);
    delete(handles.t31);delete(handles.t32);delete(handles.t41);
    delete(handles.t42);
end
figure(handles.stan);
subplot(2,2,1);
handles.s11=stem(handles.Imincs*1e3,handles.Pout(97)*1e3,'fill',':.r','LineWidth',2);
handles.t11=text(handles.Imincs*1e3,handles.Pout(97)*1e3/2,'Imin',...
                'HorizontalAlignment','right','FontSize',10);
handles.s12=stem(handles.Imaxcs*1e3,handles.Pout(97)*1e3,'fill',':.b','LineWidth',2);
handles.t12=text(handles.Imaxcs*1e3,handles.Pout(97)*1e3/2,...
                'Imax','HorizontalAlignment','left','FontSize',10);
subplot(2,2,2);
handles.s21=stem(handles.Imincs*1e3,handles.Nw(95)*1.1,'fill',':.r','LineWidth',2);
handles.t21=text(handles.Imincs*1e3,handles.Nw(95)/2,'Imin',...
                'HorizontalAlignment','right','FontSize',10);
handles.s22=stem(handles.Imaxcs*1e3,handles.Nw(95)*1.1,'fill',':.b','LineWidth',2);
handles.t22=text(handles.Imaxcs*1e3,handles.Nw(95)/2,...
                'Imax','HorizontalAlignment','left','FontSize',10);
subplot(2,2,3);
handles.s31=stem(handles.Imincs*1e3,handles.Ng(95)*1.1,'fill',':.r','LineWidth',2);
handles.t31=text(handles.Imincs*1e3,handles.Ng(95)/2,'Imin',...
                'HorizontalAlignment','right','FontSize',10);
handles.s32=stem(handles.Imaxcs*1e3,handles.Ng(95)*1.1,'fill',':.b','LineWidth',2);
handles.t32=text(handles.Imaxcs*1e3,handles.Ng(95)/2,...
                'Imax','HorizontalAlignment','left','FontSize',10);
subplot(2,2,4);
handles.s41=stem(handles.Imincs*1e3,handles.Nsch(95)*1.1,'fill',':.r','LineWidth',2);
handles.t41=text(handles.Imincs*1e3,handles.Nsch(95)/2,'Imin',...
                'HorizontalAlignment','right','FontSize',10);
handles.s42=stem(handles.Imaxcs*1e3,handles.Nsch(95)*1.1,'fill',':.b','LineWidth',2);
handles.t42=text(handles.Imaxcs*1e3,handles.Nsch(95)/2,...
                'Imax','HorizontalAlignment','left','FontSize',10);
if handles.val==1 && handles.case==1
    delete(handles.s11);delete(handles.s12);delete(handles.s21);
    delete(handles.s22);delete(handles.s31);delete(handles.s32);
    delete(handles.s41);delete(handles.s42);delete(handles.t11);
    delete(handles.t12);delete(handles.t21);delete(handles.t22);
    delete(handles.t31);delete(handles.t32);delete(handles.t41);
    delete(handles.t42);
    handles.val=0; handles.case=0;
else
    handles.val=1;
end
guidata(gcbo,handles);
%-------------------------------------------------------------------------------------------------------%
%Analisi Electro - Optic, usa il metodo di Runge-Kutta
%I grafici ottenuti riguardano I-Pout, I-Carrier Density, I-Frequency shift curves
function RunD1_Callback(hObject, eventdata, handles)

DummyVariable(eventdata);
set(get(hObject,'parent'),'Pointer','watch');
h=waitbar(1/4,'Please waiting... calculation in progress');drawnow;
%Intervallo temporale di integrazione
T=[handles.Tstart-15e-9 handles.Tend];                
%Valori iniziali del vettore di stato
vin=[1e7; 1e7; 1e7; 1e7];
waitbar(2/4);
%Modifica dei parametri di integrazione usati dalla funzione ode45
options=odeset('initialstep',1e-3, 'AbsTol',1e-6, 'RelTol', 1e-3,'MaxStep',0.1e-10);
[t,y]=ode15s(@Ele,T, vin,options,handles,1);
Nsch=y(:,1);                 % The spatially average carrier numbers in SCH
Ng=y(:,2);                   % The spatially average carrier numbers in gateway
Nw=y(:,3);                   % The spatially average carrier numbers in well
S=y(:,4);                    % The photon number in the cavity laser
waitbar(3/4);
[I]=Curre_Sig(t,handles,handles.title);     %The current signal
Pout=(handles.n0*handles.h*handles.v*S/handles.Taup);  %Potenza di uscita
waitbar(4/4);
close(h);
%Esecuzione grafici
GraficiD(t,I,Pout,Nw,Ng,Nsch,handles);
%Vettore di t, dopo (t0-0.5ns) prima del gradino
tn = handles.Tst-0.5e-9;                        % Tempo prima d'inietare il gradino
pos = t>=tn;                                    % Possizione da creare il nuovo vettore
handles.tn0 = t(pos);                           % Vettore dopo la possizione di tn
handles.Ppos = Pout(pos);                       % Vettore di potenza dopo la possizione di tn
handles.Ipos = I(pos);                          % Vettore di corrente dopo la possizione di tn
handles.Nw = Nw(pos);                   
handles.Ng = Ng(pos);                    
handles.Nsch = Nsch(pos);               
%Memorizza i dati nell'interfaccia
guidata(gcbo,handles);
set(handles.AM,'Visible','on'); %set(handles.FM,'Visible','on');
set(handles.AM,'Position',[10,12.308,23.8,2]);
%set(handles.FM,'Position',[10,9.385,23.8,2]);
set(get(hObject,'parent'),'Pointer','arrow');
%-------------------------------------------------------------------------------------------------------%
%Analisi Switching On
function RunD2_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
set(get(hObject,'parent'),'Pointer','watch');
h=waitbar(1/4,'Please waiting... calculation in progress');drawnow;
%Intervallo temporale di integrazione
T=[handles.Tstart-2e-9 handles.Tend];
%Valori iniziali del vettore di stato
vin=[1; 1; 1; 1];
waitbar(2/4);
%Modifica dei parametri di integrazione usati dalla funzione ode45
options=odeset('initialstep',1e-3, 'AbsTol',1e-6, 'RelTol', 1e-3,'MaxStep',0.1e-10);
[t,y]=ode15s(@Ele,T, vin,options,handles,handles.title);
Nsch=y(:,1);                 % The spatially average carrier numbers in SCH
Ng=y(:,2);                   % The spatially average carrier numbers in gateway
Nw=y(:,3);                   % The spatially average carrier numbers in well
S=y(:,4);                    % The photon number in the cavity laser
waitbar(3/4);
[I]=Curre_Sig(t,handles,handles.title);     %The current signal
Pout=(handles.n0*handles.h*handles.v*S/handles.Taup);  %Potenza di uscita
waitbar(4/4);
close(h);
%Esecuzione grafici
GraficiD(t,I,Pout,Nw,Ng,Nsch,handles);
%Vettore di t, dopo (t0-0.5ns) prima del gradino
tn = handles.Tst-0.5e-9;                        % Tempo prima d'inietare il gradino
pos = t>=tn;                                    % Possizione da creare il nuovo vettore
handles.tn0 = t(pos);                           % Vettore dopo la possizione di tn
handles.Ppos = Pout(pos);                       % Vettore di potenza dopo la possizione di tn
handles.Ipos = I(pos);                          % Vettore di corrente dopo la possizione di tn
handles.Nw = Nw(pos);                   
handles.Ng = Ng(pos);                    
handles.Nsch = Nsch(pos);               
%Memorizza i dati nell'interfaccia
guidata(gcbo,handles);
set(get(hObject,'parent'),'Pointer','arrow');
%-------------------------------------------------------------------------------------------------------%
%Analisi in AM
function AM_Callback(hObject, eventdata, handles)
%Interpolazione a passo costante, prima tn0=0 ns.
analisi='AM';
DummyVariable(hObject,eventdata);
handles.TStart = (handles.tn0-handles.tn0(1))*1e9;                      % Tempo d'inizio, in 0 in ns
handles.TimeStep = 1e-4;                                                % Tempo di campionamento
handles.Tint = handles.TStart(1):handles.TimeStep:handles.TStart(end);  % Tempo con paso costante
handles.Poutint = interp1(handles.TStart,handles.Ppos,handles.Tint);    % Potenza interpolata costante
handles.Iint = interp1(handles.TStart,handles.Ipos,handles.Tint);       % Corrente interpolata costante
Campo=sqrt(handles.Poutint);
%Chiamata alla funzione del gradino
IRGradino(handles.Tint,handles.Iint,Campo,analisi);
%-------------------------------------------------------------------------------------------------------%
%Analisi in FM
function FM_Callback(hObject, eventdata, handles)
analisi='FM';
DummyVariable(hObject,eventdata);
handles.TStart = (handles.tn0-handles.tn0(1))*1e9;                      % Tempo d'inizio, in 0 in ns
handles.TimeStep = 1e-4;                                                % Tempo di campionamento
handles.Tint = handles.TStart(1):handles.TimeStep:handles.TStart(end);  % Tempo con paso costante
handles.Stint = interp1(handles.TStart,handles.S,handles.Tint);  % Potenza interpolata costante
handles.Iint = interp1(handles.TStart,handles.Ipos,handles.Tint);       % Corrente interpolata costante
Campo=sqrt((handles.Stint)*1e-9);
%Chiamata alla funzione del gradino
IRGradino(handles.Tint,handles.Iint,Campo,analisi);
%-------------------------------------------------------------------------------------------------------%
function param=IRSetup(TSimStart,TSimStop,TSimStep)
param.Start=[0 1];
TStartSample=ceil((TSimStart)/TSimStep)+1;
TStopSample=floor((TSimStop)/TSimStep)+1;
ptTot=TStopSample-TStartSample+1;
deltaPt=2^floor(log2(ptTot*0.05));
ptOpt=floor(ptTot/deltaPt)*deltaPt;
param.Points=ptOpt;
param.Stop=[ptOpt*TSimStep,(ptOpt+1),TSimStop];
%-------------------------------------------------------------------------------------------------------%
%Calcolo della risposta al gradino 
function IRGradino(Tem,Curr,Campo1,AN)

%Impostazione variabili
TStart=0;
TEnd=Tem(end);
TimeStep=Tem(2)-Tem(1);
%Funzione di attesa
h=waitbar(1/6,'Please waiting... calculation in progress');drawnow;
IRP=IRSetup(TStart,TEnd,TimeStep);
if Curr(IRP.Start(2))~=Curr(IRP.Start(2)+1) &&  Curr(IRP.Start(2)+1)==Curr(IRP.Stop(2))
   IRP.Start(2)=IRP.Start(2)+1;
end
campo=Campo1(IRP.Start(2):IRP.Stop(2));
PositiveIndex=ceil(IRP.Points/2+1):IRP.Points;
FreqRange=linspace(-0.5/TimeStep,0.5/TimeStep,IRP.Points);
FreqRange=FreqRange(PositiveIndex);
waitbar(2/6);
FFT.Potenza=fftshift(fft(gradient(abs(campo).^2,TimeStep)))/IRP.Points*(IRP.Stop(1)-IRP.Start(1));
FFT.Potenza=FFT.Potenza(PositiveIndex);
waitbar(3/6);
if AN=='AM'
    IR.AM=FFT.Potenza/FFT.Potenza(2);
    waitbar(4/6);
    clear campo;
    clear corrente;
    waitbar(5/6);
    %Grafico
    figure;
    semilogx(FreqRange,20*log10(abs(IR.AM)));
    waitbar(6/6);
    close (h);
    ylabel('Normalized AM Impulse Response [dB]');
    xlabel('Frequency [GHz]');
    grid on;
elseif AN=='FM'
    IR.FM_NN=FFT.Potenza;
    IR.FM=FFT.Potenza/FFT.Potenza(4);
    waitbar(4/6);
    clear campo;
    clear corrente;
    waitbar(5/6);
    %Grafici
    figure;
    semilogx(FreqRange,20*log10(abs(IR.FM)));
    waitbar(6/6);
    close (h);
    ylabel('Normalized FM Impulse Response [dB]');
    xlabel('Frequency [GHz]');
    grid on;
end
%-------------------------------------------------------------------------------------------------------%
%Analisi PRBS, usa il metodo Runge-Kutta
function RunD3_Callback(hObject, eventdata, handles)

DummyVariable(hObject,eventdata);
set(get(hObject,'parent'),'Pointer','watch');
h=waitbar(1/6,'Please waiting... calculation in progress');drawnow;
%Intervallo temporale di integrazione
T=[handles.Tstart-15e-9 (handles.Nbit/handles.Rb)+15e-9];        
%Valori iniziali del vettore di stato
vin=[1e10; 1e10; 1e10; 1e10];
waitbar(2/4);
%Modifica dei parametri di integrazione usati dalla funzione ode45
options=odeset('initialstep',1e-3,'AbsTol',1e-6,'RelTol',1e-3,'MaxStep',1e-11);
[t,y]=ode15s(@Ele,T, vin,options,handles,handles.title);
Nsch=y(:,1);                 % The spatially average carrier numbers in SCH
Ng=y(:,2);                   % The spatially average carrier numbers in gateway
Nw=y(:,3);                   % The spatially average carrier numbers in well
S=y(:,4);                    % The photon number in the cavity laser
waitbar(3/4);
[I]=Curre_Sig(t,handles,handles.title);     %The current signal
Pout=(handles.n0*handles.h*handles.v*S/handles.Taup);  %Potenza di uscita
CER = (handles.Imax-1)/(handles.Imin-1);
handles.CER_dB = 10*log10(abs(CER));
waitbar(4/4);
close(h);
%Esecuzione grafici
GraficiD(t,I,Pout,Nw,Ng,Nsch,handles);
%Vettore di t, dopo (t0-0.5ns) prima del gradino
tn = handles.Tst-0.5e-9;                        % Tempo prima d'inietare il gradino
pos = t>=tn;                                    % Possizione da creare il nuovo vettore
handles.tn0 = t(pos);                           % Vettore dopo la possizione di tn
handles.Ppos = Pout(pos);                       % Vettore di potenza dopo la possizione di tn
handles.Ipos = I(pos);                          % Vettore di corrente dopo la possizione di tn
handles.Nw = Nw(pos);                   
handles.Ng = Ng(pos);                    
handles.Nsch = Nsch(pos);               
%Memorizza i dati nell'interfaccia
guidata(gcbo,handles);
set(handles.P_CER,'Visible','on'); set(handles.T_CER,'Visible','on');
set(handles.R_CER,'Visible','on'); set(handles.R_CER,'String',num2str(handles.CER_dB,'%.2f'));
set(get(hObject,'parent'),'Pointer','arrow');
%-------------------------------------------------------------------------------------------------------%
%Equazione differenziale per analisi segnale
function out=Ele(t,invet,handles,tipo)
%Richiamo la funzione P_I
[I]=Curre_Sig(t,handles,tipo);
%Per semplificare le espressioni, metto nella variabile N ed Np i valori del vettore di stato
handles.Nsch = invet(1);            % The spatially averaged carrier number in sch
handles.Ng = invet(2);              % The spatially averaged carrier number in gateway
handles.Nw = invet(3);              % The spatially averaged carrier number in well
handles.S = invet(4);               % The photon number in the laser cavity
handles.Nwn=handles.Nw/handles.Vw;
handles.go = handles.a0*(handles.Nwn-handles.Nwo);
%Rate equations che definiscono
handles.dNsch=(handles.ni*I/handles.q)-(handles.Nsch*(1/handles.Tauns+1/...
               handles.Taud))+(handles.Ng/handles.Taug);
handles.dNg=(handles.Nsch/handles.Taud)-(handles.Ng*(1/handles.Taucap+1/...
             handles.Taug+1/handles.Taung))+(handles.Nw/handles.Tauesc);
handles.dNw=(handles.Ng/handles.Taucap)-(handles.Nw*(1/handles.Tauesc+1/...
             handles.Taunw))-(handles.Vg*handles.go*handles.Gamw*...
            (handles.S/(1+handles.Ep*handles.Gamw*handles.S+handles.Ep2*...
            (handles.Gamw*handles.S).^2)));
handles.dS=handles.Espont+handles.S*(handles.Gamw*handles.Vg*handles.go./...
          (1+handles.Ep*handles.Gamw*handles.S+handles.Ep2*(handles.Gamw*...
           handles.S).^2)-1/handles.Taup);
out=[handles.dNsch;handles.dNg;handles.dNw;handles.dS];
%Memorizza i dati nell'interfaccia
guidata(gcbo,handles);
%-------------------------------------------------------------------------------------------------------%
%Simula il comportamento del segnale di corrente iniettata nel laser con il rispetto del tempo
function [I]=Curre_Sig(t,handles,tipo1)
switch(tipo1)
    case 1
        %Step di corrente
        I = handles.Imincs*(t>=handles.Tstart-15e-9)+handles.Imax*(t>=handles.Tst);            
        guidata(gcbo,handles);
    case 2
        %Gradino
        I = handles.Imincs*(t>=handles.Tstart-2e-9)+(handles.Imaxcs-handles.Imincs)*(t>=handles.Tst);  
        guidata(gcbo,handles);  
    case 3
        t=t+15e-9;
        Tb = 1/handles.Rb; % Tempo di bit
        Ipico = handles.Imaxcs-handles.Imincs;                % Ampiezza picco picco della modulazione
        vett_PRBS = prbs(handles.Nbit);
        vett_PRBS = [ones(1,(11e-9*handles.Rb)),0,0,vett_PRBS];
        I=zeros(1,length(t));
        for k=1:length(t)
            I(k) = handles.Imincs*(t(k)>=handles.Tstart)+pulse_PRBS(t(k),vett_PRBS,Tb,4,Ipico)*(t(k)>=handles.Tst);
        end
        guidata(gcbo,handles);  
end
%-------------------------------------------------------------------------------------------------------%
%Funzione usata per la generazione dei parametri di corrente
function [Ith]=GetSoglia(I,P)
warning('off','MATLAB:polyfit:RepeatedPoints');
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
PoNmax=max(P);
Validi=(P>PoNmax/100);
P=P(Validi);
I=I(Validi);
IndiceMax=min([30,length(P)]);
P=P(1:IndiceMax);
I=I(1:IndiceMax);
c=polyfit(P,I,3);
Ith=c(end);
%-------------------------------------------------------------------------------------------------------%
function z=prbs(N)
s = rand('state');
rand('state',0);
init=rand(1,N);
rand('state',s);
init(init>0.5)=1;
init(init~=1)=0;
g=[N 1];
n=length(init);
z=zeros(1,2^n-1);
z(1:n)=init;
for i=(n+1):(2^n-1)
    q=z(i-g(1));
    for j=2:length(g)
        q=xor(q,z(i-g(j)));
    end
    z(i)=q;
end
%-------------------------------------------------------------------------------------------------------%
function  pbit=pulse_PRBS(t,seq,T,x,Peak)

ib=ceil(t/T);
if ib==0
   pbit=0;
   return
else
    if (ib~=1)
        if (ib<=length(seq))
            if (seq(ib)==1) && (seq(ib-1)==0)
                tval=t-(ib-1)*T-T/2*(1+1/x);
                p=fun_rcos(tval,1/T,Peak,1/x);
            elseif (seq(ib)==1) && (seq(ib-1)==1)
                p=Peak;
            elseif (seq(ib)==0) && (seq(ib-1)==1)
                tval=t-(ib-1)*T-T/2*(1+1/x)+T;
                p=fun_rcos(tval,1/T,Peak,1/x);
            elseif (seq(ib)==0) && (seq(ib-1)==0)
                p=0;
            end
        else
            p=0;
        end
    else   % primo bit
        if seq(ib)==0
            p=0;
        else
            tval=t-(ib-1)*T-T/2*(1+1/x);
            p=fun_rcos(tval,1/T,Peak,1/x);
        end
    end
    pbit=p;
end
%-------------------------------------------------------------------------------------------------------%
function  fc=fun_rcos(fin,Ts,Peak,alpha_r)

fc=zeros(1,length(fin));
for ifr=1:length(fin)
   f=fin(ifr);
   if (abs(f)<((1-alpha_r)/(2*Ts)))
       fc(ifr)=Peak;
   elseif ((abs(f)>=((1-alpha_r)/(2*Ts))) && (abs(f)<=((1+alpha_r)/(2*Ts))))
       fc(ifr)=Peak/2*(1- sin(pi*Ts/alpha_r*(abs(f)-1/(2*Ts))));
   else
       fc(ifr)=0;
   end
end
%-------------------------------------------------------------------------------------------------------%
function GraficiA(i,handles)
leg={};
subplot(2,2,1);
plot(handles.I*1e3,handles.Poutn*1e3,'LineWidth',1.5);
hold on;
leg{i}=[handles.quantitascansione ' = ',num2str(handles.ValorePerlegenda)];
legend(leg,4);
grid on
xlabel('Current [mA]');
ylabel('Output Power [mW]');
title('P vs I');
gradeplots;
subplot(2,2,2);
plot(handles.I*1e3,handles.Nw,'LineWidth',1.5);
hold on;
leg{i}=[handles.quantitascansione ' = ',num2str(handles.ValorePerlegenda)];
legend(leg,4);
grid on
xlabel('Current [mA]');
ylabel('Carrier Density in quantum well [cm-3]');
title('Carrier Density in quantum Well vs I');
gradeplots;
subplot(2,2,3);
plot(handles.I*1e3,handles.Ng,'LineWidth',1.5);
hold on;
leg{i}=[handles.quantitascansione ' = ',num2str(handles.ValorePerlegenda)];
legend(leg,4);
grid on
xlabel('Current [mA]');
ylabel('Carrier density in the gateway [cm-3]');
title('Carrier density in gateway vs I');
gradeplots;
subplot(2,2,4);
plot(handles.I*1e3,handles.Nsch,'LineWidth',1.5);
hold on;
leg{i}=[handles.quantitascansione ' = ',num2str(handles.ValorePerlegenda)];
legend(leg,4);
grid on
xlabel('Current [mA]');
ylabel('Carrier density in the SCH [cm-3]');
title('Carrier density in SCH vs I');
gradeplots;
guidata(gcbo,handles);
%-------------------------------------------------------------------------------------------------------%
function GraficiD(Gt,GI,GPout,GNw,GNg,GNsch,handles)

set(0,'Units','pixels'); 
scnsize = get(0,'ScreenSize');
pos=[5,35,scnsize(3)-8,scnsize(4)-115];
switch (handles.title)
    case 1
        figure('Name','Dynamic analysis - Electro-Optic','NumberTitle','off','Position',pos);
        tn0 = handles.Tst-1e-9;
        Gxmin=handles.Tst*1e9-1;
        Gxmax=handles.Tend*1e9;
    case 2
        figure('Name','Dynamic analysis - Switching On','NumberTitle','off','Position',pos);       
        tn0 = handles.Tst-1e-9;
        Gxmin=handles.Tst*1e9-1;
        Gxmax=handles.Tend*1e9;
    case 3
        figure('Name','Dynamic analysis - PRBS Signal','NumberTitle','off','Position',pos);
        tn0 = handles.Tst;
        Gxmin=handles.Tst*1e9;
        Gxmax=((handles.Nbit/handles.Rb)+handles.Tst)*1e9;
end
p = Gt>=tn0;        % Possizione da creare il nuovo vettore
Gtp1 = Gt(p);       % Vettore dopo la possizione di tn
GIp1 = GI(p);       % Vettore di corrente dopo la possizione di tn
GPp1 = GPout(p);    % Vettore di potenza dopo la possizione di tn
GNwp1 = GNw(p);
GNgp1 = GNg(p);
GNschp1 = GNsch(p);
subplot(2,2,1,'v6');
[ax,h1,h2] = plotyy(Gtp1*1e9,GIp1*1e3,Gtp1*1e9,GPp1*1e3,'plot');
grid on;
set(get(ax(1),'Ylabel'),'String','Current [mA]');
set(get(ax(2),'Ylabel'),'String','Output Power [mw]');
xlabel('Time [ns]');
set(ax,'xlim',[Gxmin Gxmax]);
title('P_{out} - Current vs Time');
set(h1,'LineStyle',':','color','b','LineWidth',2.5);
set(h2,'LineStyle','-','color','r','LineWidth',2.5);
set(ax(1),'YAxisLocation','r','Ycol','b');
set(ax(2),'YAxisLocation','l','Ycol','r','YMinorGrid','on');
subplot(2,2,2,'v6');
[ax,h1,h2] = plotyy(Gtp1*1e9,GIp1*1e3,Gtp1*1e9,GNwp1,'plot');
grid on;
set(get(ax(1),'Ylabel'),'String','Current [mA]');
set(get(ax(2),'Ylabel'),'String','Carrier Density in quantum well [cm^-3]');
xlabel('Time [ns]');
set(ax,'xlim',[Gxmin Gxmax]);
title('Carrier Density in quantum well - Currente vs Time');
set(h1,'LineStyle',':','color','b','LineWidth',2.5);
set(h2,'LineStyle','-','color','r','LineWidth',2.5);
set(ax(1),'YAxisLocation','r','Ycol','b');
set(ax(2),'YAxisLocation','l','Ycol','r','YMinorGrid','on');
subplot(2,2,3,'v6');
[ax,h1,h2] = plotyy(Gtp1*1e9,GIp1*1e3,Gtp1*1e9,GNgp1,'plot');
grid on;
set(get(ax(1),'Ylabel'),'String','Current [mA]');
set(get(ax(2),'Ylabel'),'String','Carrier Density in the gateway [cm^-3]');
xlabel('Time [ns]');
set(ax,'xlim',[Gxmin Gxmax]);
title('Carrier Density in the gateway - Current vs Time');
set(h1,'LineStyle',':','color','b','LineWidth',2.5);
set(h2,'LineStyle','-','color','r','LineWidth',2.5);
set(ax(1),'YAxisLocation','r','Ycol','b');
set(ax(2),'YAxisLocation','l','Ycol','r','YMinorGrid','on');
subplot(2,2,4,'v6');
[ax,h1,h2] = plotyy(Gtp1*1e9,GIp1*1e3,Gtp1*1e9,GNschp1,'plot');
grid on;
set(get(ax(1),'Ylabel'),'String','Current [mA]');
set(get(ax(2),'Ylabel'),'String','Carrier Density in the SCH [cm^-3]');
xlabel('Time [ns]');
set(ax,'xlim',[Gxmin Gxmax]);
title('Carrier Density in the SCH - Currente vs Time');
set(h1,'LineStyle',':','color','b','LineWidth',2.5);
set(h2,'LineStyle','-','color','r','LineWidth',2.5);
set(ax(1),'YAxisLocation','r','Ycol','b');
set(ax(2),'YAxisLocation','l','Ycol','r','YMinorGrid','on');