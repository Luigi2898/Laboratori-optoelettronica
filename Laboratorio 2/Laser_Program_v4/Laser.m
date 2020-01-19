%-------------------------------------------------------------------------%
function varargout = Laser(varargin)
% LASER M-file for Laser.fig
%Funzione per la gestione dell'interfaccia di laser
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Laser_OpeningFcn, ...
                   'gui_OutputFcn',  @Laser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code
%-------------------------------------------------------------------------%
function Laser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% Choose default command line output for Laser
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
%-------------------------------------------------------------------------%
function varargout = Laser_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;
%-------------------------------------------------------------------------%
function MENU_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
%-------------------------------------------------------------------------%
% It is used to suppress warning from M-Lint for variable not
% used and functions never called. The formers are caused by parameters
% automatically passed by MatLab in callback functions; the latters are
% functions used as callback and therefore passed as a string and
% automatically called by MatLab
function DummyVariable(varargin)
if(0)
   MENU_Callback;
   OK_Callback;
end
%-------------------------------------------------------------------------%
function MENU_Callback(hObject, eventdata, handles)
DummyVariable(hObject,eventdata);
%struttura case per lo smistamento a seconda della scelta iniziale del modello di analisi
switch get(handles.MENU,'Value')
    case 1
        set(handles.T2,'String','');
        cla;
    case 2
        %scelta del modello a 1 portatore                
        set(handles.T2,'String',{'The bulk semiconductor laser model is mathematically described with two differential equation:',...
                                 '','     > dN/dt - carriers density',...
                                 '','     > dNp/dt - photon density'});
        %Visualizzazione immagine 
        axes(handles.Image);
        image(imread('Population1.jpg'));
        axis equal;
        axis tight;
        axis off;
    case 3
        %scelta del modello a 2 portatori
        set(handles.T2,'String',{'The model of the quantum well material is mathematically described with three differential equation:',...
                                 '','     > dNb/dt - carriers density in the barrier',...
                                 '','     > dN/dt - carriers density in the well',...
                                 '','     > dNp/dt - photon density'});
        %Visualizzazione immagine 
        axes(handles.Image);
        image(imread('Population2.jpg'));
        axis equal;
        axis tight;
        axis off;
    case 4
        %scelta del modello a 3 portatori
        set(handles.T2,'String',{'The quantum well semiconductor laser model can be furthermore completed,assuming three populations and therefore is composed of four differential equation:',...
                                 '','     > dNsch/dt - carrier density in the separate                                  confinement heterostructure',...
                                 '','     > dNg/dt - carrier density in the gateway                                 (region between the QW)',...
                                 '','     > dNw/dt - carrier density in the well',...
                                 '','     > dS/dt - photon number in the laser cavity'});
        %Visualizzazione immagine 
        axes(handles.Image);
        image(imread('Population3.jpg'));
        axis equal;
        axis tight;
        axis off;
end
%-------------------------------------------------------------------------%
function OK_Callback(hObject, eventdata, handles)
DummyVariable(hObject,eventdata);
%Impostazione della variabile relativa al tasto di ok       
switch get(handles.MENU,'Value')
    
    case 2
    Laser1;
    
    case 3
    Laser2;
    
    case 4
    Laser3;
    
end
%-------------------------------------------------------------------------%
function CLOSE_Callback(hObject, eventdata, handles)
DummyVariable(hObject,eventdata);

clear all;
close all;
