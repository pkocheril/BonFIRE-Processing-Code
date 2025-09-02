% Dichroic/bandpass name lookup
function [NAME] = filterlookup(INPUT)
%%% This function looks up the dichroic/bandpass filter name for a specified
% numeric input (generated from LabVIEW).

    % In case nothing found
    NAME = string(INPUT);
    NAME = char(NAME);

    % Dichroics
    if INPUT == 1
        NAME = 'Ag';
    end
    if INPUT == 90
        NAME = 'BSX10R';
    end
    if INPUT == 91
        NAME = 'UFBS9010';
    end
    if INPUT == 425
        NAME = 'DMSP425R';
    end
    if INPUT == 505
        NAME = '505-SDi01';
    end
    if INPUT == 509
        NAME = '509-Di01';
    end
    if INPUT == 526
        NAME = '526-Di01';
    end
    if INPUT == 535
        NAME = '535-SDi01';
    end
    if INPUT == 594
        NAME = 'Di03-594';
    end
    if INPUT == 660
        NAME = 'Di03-R660';
    end
    if INPUT == 700
        NAME = '700-Di01';
    end
    if INPUT == 738
        NAME = '738-Di01';
    end
    if INPUT == 765
        NAME = '765-Di01';
    end
    if INPUT == 801
        NAME = '801-Di02';
    end
    if INPUT == 830
        NAME = 'Di02-R830';
    end

    % Bandpasses
    if INPUT == 439
        NAME = '439/154';
    end
    if INPUT == 454
        NAME = 'FBH450/40';
    end
    if INPUT == 475
        NAME = '475/50';
    end
    if INPUT == 550
        NAME = '550/49';
    end
    if INPUT == 585
        NAME = '585/29';
    end
    if INPUT == 609
        NAME = '609/54';
    end
    if INPUT == 610
        NAME = '609/181';
    end
    if INPUT == 650
        NAME = '650/60';
    end
    if INPUT == 665
        NAME = '665/150';
    end
    if INPUT == 681
        NAME = '681/24';
    end
    if INPUT == 708
        NAME = '708/75';
    end
    if INPUT == 709
        NAME = '709/167';
    end
    if INPUT == 719
        NAME = '719/60';
    end
    if INPUT == 731
        NAME = '731/137';
    end
    if INPUT == 769
        NAME = '769/41';
    end
    if INPUT == 775
        NAME = '775/140';
    end
    if INPUT == 792
        NAME = '792/64';
    end
    if INPUT == 893
        NAME = '893/209';
    end
    
    % Notch filters
    if INPUT == 514
        NAME = 'NF514/17';
    end
    if INPUT == 658
        NAME = 'NF658/26';
    end

    % Longpasses
    if INPUT == 400
        NAME = 'FELH0400';
    end
    
    % Shortpasses
    if INPUT == 45
        NAME = 'FESH0450';
    end
    if INPUT == 50
        NAME = 'FESH0500';
    end
    if INPUT == 60
        NAME = 'FESH0600';
    end
    if INPUT == 70
        NAME = 'FESH0700';
    end
    if INPUT == 80
        NAME = 'FESH0800';
    end
    if INPUT == 85
        NAME = 'FESH0850';
    end

    return
end
