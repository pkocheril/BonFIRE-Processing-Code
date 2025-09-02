% Calculate solvent reaction fields
function [SRF,ONSAGERFACTOR] = onsager(DIELEC,DIPOLEDEBYE,REFRAC,VOLUMEA3)
% This function calculates Onsager solvent reaction fields given a list of
% dielectric constants and the solute properties of interest (dipole moment
% in D, refractive index, and volume in A^3)

    % Set up vectors (possibly being fed a list of DIELEC)
    SRF = zeros(length(DIELEC),1); numdielec = SRF;
    ONSAGERFACTOR = SRF;

    % Parse inputs - replace empty values with Rh800
    if isempty(DIPOLEDEBYE)
        DIPOLEDEBYE = 5.3383; % solute dipole moment, D
    end
    if isempty(REFRAC)
        REFRAC = 1.5; % solute refractive index
    end
    if isempty(VOLUMEA3)
        VOLUMEA3 = 164; % solute molecular volume, A^3
    end

    % Set up constants
    c = 299792458; % m/s
    permit = 8.85418782e-12; % F/m, permittivity of free space
    coulomb = 1/(4*pi*permit); % N m^2 C^-2, Coulomb constant
    % h = 6.62607015e-34; % J s or kg m^2 s^-1
    % hbar = h/(2*pi); % J s
    % kB = 1.380649e-23; % J/K
    % temperature = 300; % K
    % kT = kB*temperature; % J
    % c100 = c*100; % cm/s, or Hz per cm-1
    % NAv = 6.02214076e23; % molecules per mole
    % amukg = 1/(NAv*1e3); % kg per 1 amu (bc 1 mol of 1 amu particles weighs 1 g)
    % angm = 1e-10; % m per 1 angstrom

    % Calculate SRFs
    dipoleSI = DIPOLEDEBYE*1e-21/c; % dipole moment in C m
    volumeSI = VOLUMEA3*1e-30; % molecular volume in m^3

    % Convert to numbers if needed
    for i=1:length(DIELEC)
        if isstring(DIELEC(i)) || ischar(DIELEC(i)) % look up dielectric constants if fed solvent names
            numdielec(i) = 1; % obviously wrong if unmatched

            % Organics
            if strcmpi(DIELEC(i),"DMF") || strcmpi(DIELEC(i),"dimethylformamide")
                numdielec(i) = 36.7;
            end
            if strcmpi(DIELEC(i),"DOX") || strcmpi(DIELEC(i),"dioxane") || strcmpi(DIELEC(i),"1,4-dioxane")
                numdielec(i) = 2.25;
            end
            if strcmpi(DIELEC(i),"EtOAc") || strcmpi(DIELEC(i),"ethyl acetate")
                numdielec(i) = 6.02;
            end
            if strcmpi(DIELEC(i),"THF") || strcmpi(DIELEC(i),"tetrahydrofuran")
                numdielec(i) = 7.52;
            end
            if strcmpi(DIELEC(i),"Tol") || strcmpi(DIELEC(i),"toluene")
                numdielec(i) = 2.38;
            end
            if strcmpi(DIELEC(i),"Hex70") || strcmpi(DIELEC(i),"7:3 hexane:chloroform")
                numdielec(i) = 3.09;
            end
            if strcmpi(DIELEC(i),"CHCl3") || strcmpi(DIELEC(i),"chloroform")
                numdielec(i) = 4.81;
            end
            if strcmpi(DIELEC(i),"BBZ") || strcmpi(DIELEC(i),"benzyl benzoate")
                numdielec(i) = 4.9;
            end
            if strcmpi(DIELEC(i),"DBE") || strcmpi(DIELEC(i),"dibenzyl ether")
                numdielec(i) = 3.86;
            end

            % Water/DMSO mixtures
            if strcmpi(DIELEC(i),"PBS") || strcmpi(DIELEC(i),"PBS100") || strcmpi(DIELEC(i),"water") || strcmpi(DIELEC(i),"D2O") || strcmpi(DIELEC(i),"H2O")
                numdielec(i) = 80.1;
            end
            if strcmpi(DIELEC(i),"PBS80") || strcmpi(DIELEC(i),"DMSO20") || strcmpi(DIELEC(i),"D2O80")
                numdielec(i) = (0.8*80.1)+(0.2*46.7);
            end
            if strcmpi(DIELEC(i),"PBS60") || strcmpi(DIELEC(i),"DMSO40") || strcmpi(DIELEC(i),"D2O60")
                numdielec(i) = (0.6*80.1)+(0.4*46.7);
            end
            if strcmpi(DIELEC(i),"PBS40") || strcmpi(DIELEC(i),"DMSO60") || strcmpi(DIELEC(i),"D2O40")
                numdielec(i) = (0.4*80.1)+(0.6*46.7);
            end
            if strcmpi(DIELEC(i),"PBS20") || strcmpi(DIELEC(i),"DMSO80") || strcmpi(DIELEC(i),"D2O20")
                numdielec(i) = (0.2*80.1)+(0.8*46.7);
            end
            if strcmpi(DIELEC(i),"PBS0") || strcmpi(DIELEC(i),"DMSO") || strcmpi(DIELEC(i),"DMSO100") || strcmpi(DIELEC(i),"dimethylsulfoxide")
                numdielec(i) = 46.7;
            end
            
            % Alcohols
            if strcmpi(DIELEC(i),"MeOH") || strcmpi(DIELEC(i),"methanol")
                numdielec(i) = 32.7;
            end
            if strcmpi(DIELEC(i),"EtOH") || strcmpi(DIELEC(i),"ethanol")
                numdielec(i) = 24.5;
            end
            if strcmpi(DIELEC(i),"iPrOH") || strcmpi(DIELEC(i),"isopropanol") || strcmpi(DIELEC(i),"isopropyl alcohol")
                numdielec(i) = 17.9;
            end
            if strcmpi(DIELEC(i),"BuOH") || strcmpi(DIELEC(i),"butanol") || strcmpi(DIELEC(i),"n-butanol")
                numdielec(i) = 17.8;
            end
            if strcmpi(DIELEC(i),"BnOH") || strcmpi(DIELEC(i),"BzOH") || strcmpi(DIELEC(i),"benzyl alcohol")
                numdielec(i) = 13;
            end
        else % given numerical value
            numdielec(i) = DIELEC(i);
        end
    end

    for i=1:length(numdielec)
        srfpart1 = -1*dipoleSI/(volumeSI); % 
        srfpart2 = (2/3).*(numdielec(i)-1).*(REFRAC^2+2)./(2.*numdielec(i)+REFRAC^2);
        SRF(i) = srfpart1.*srfpart2*coulomb*1e-8; % MV/cm
        ONSAGERFACTOR(i) = -1*srfpart2; % without molecule-specific part
    end
    return
end
