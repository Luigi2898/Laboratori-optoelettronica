% load_data: caricamento dei dati di ingresso 

[filegeo,pathgeo]=uigetfile('*.in','Choose file with the data to be loaded ');
DATA=textread([pathgeo filegeo],'%s','commentstyle','matlab','whitespace','\n');
for idata=1:length(DATA)
   eval(DATA{idata});
end