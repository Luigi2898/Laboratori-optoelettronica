% Simulazione comportamento statico SOA in MQW:  PARTE I: proprietà del
% materiale
% M. Gioannini - Maggio 2007



global material
global const_phys
global guess_Fermi


define_constants

% -------------------------------------------------------------------------
% scelta del materiale del SOA

% caricamento del file con parametri del materiale SOA

% load_data: caricamento dei dati di ingresso 

[filegeo,pathgeo]=uigetfile('*.in','Choose file with material parameters');
DATA=textread([pathgeo filegeo],'%s','commentstyle','matlab','whitespace','\n');
for idata=1:length(DATA)
   eval(DATA{idata});
end

% -------------------------------------------------------------------------

QW_properties
% calcolo di proprietà QW

calc_material_optical_prop