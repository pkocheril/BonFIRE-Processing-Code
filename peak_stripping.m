%% Part 1 of Raman baseline correction
function [Baseline, stripping]=peak_stripping(Spectrum,Window)
    stripping=0;
    y=sgolayfilt(Spectrum,0,Window);
    Baseline=zeros(length(Spectrum),1);
    for i=1:1:length(Spectrum)
       if Spectrum(i)>y(i)
           stripping=1;
           Baseline(i)=y(i);
       else
           Baseline(i)=Spectrum(i);
       end  
    end
    return
end