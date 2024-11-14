function [ GrainData,cs_all, PhaseData] = fGTranslate(GrainData, PhaseData)
%fGTranslate - translates the names of phases from astro to MTEX to load
%the correct cif file.
%CMM 2019
%Utilises MTEX codes - working with v5
for phasen = 1:GrainData.num_phase
   
    %Here is a library that compares names from BRUKER to MTEX 5.1.1 and up
    %There are some phases which have problems: fcc Moly, and hcp Nickel
    %These are presently not included
    
    %Here's how to add phases: 
    %put in the name from bruker/oxford with quotes, and a comma, and the name of the cif file
    phasenames={'Silver' , 'Silver'; 
    'Corundum' , 'Corundum';
    'Aluminium' , 'Aluminium';
    'Aluminum' , 'Aluminium';
    'Diamond', 'Diamond';
    'Chromium' , 'Chromium';
    'Chromium - alpha' ,  'Chromium';
    'Copper', 'Copper';
    'Iron', 'Iron-alpha';
    'Austenite, fcc (New)', 'Iron-Austenite';
    'Ferrite, bcc (New)', 'Iron-alpha';
    'Iron-beta', 'Iron-beta';
    'Iron-delta', 'Iron-delta';
    'Hematite', 'Hematite';
    'Magnetite', 'Magnetite';
    'Magnesium' , 'Magnesium';
    'Manganese-gamma','Manganese-gamma';
    'Molybdenum','Molybdnenum-bcc'; %(CATCH IF FCC)
    'Halite','Halite';
    'GaAs','GaAs';
    'Nickel','Nickel-cubic';
    %'Nickel','Nickel-hexagonal'; %(CATCH Is this a thing??)
    'Olivine','olivin';
    'Quartz','quartz';
    'Silicon','Silicon';
    'Titanium','Titanium-alpha';
    'Titanium-beta','Titanium-beta';
    'Baddeleyite','Baddeleyite'
    'Zirconium','Zirconium-bcc';
    'Zirconium - alpha','Zirconium-hcp';
    'Moissanite', 'Moissanite';
    'Tungsten', 'Tungsten';
    'Lithium','Lithium';
    'Forsterite','Forsterite';
    'Cementit (New)','Cementit (New)';
    'Zirc Hydride (FCC)','ZircHydride(delta)'};% PK 20200309
    
 
    phasenamenumber = strcmp([phasenames(:,1)], PhaseData(phasen).Name); % find the row with the right name
    GrainData.cifs{phasen}=char(phasenames(phasenamenumber,2)); %use the other column to translate
    
    CIF_name=fullfile(PhaseData.CIFLoc,GrainData.cifs{phasen});
    
    cs=loadCIF(CIF_name); %load from the AstroCIFs

    PhaseData(phasen).Name=char(phasenames(phasenamenumber,2));%Also translate PhaseData names for consistency
    
    if sum(phasenamenumber)>1 %Check if the names are indistinguishable - this will break lBurgLine in particular.
        error('The phase names are the same - please fix this in the Bruker file')
    end
    
    if phasen == 1
        cs_all=cs;
    else
        cs_all={cs_all,cs}; %put it all into a variable
    end
end


end
