% Pairwise statistical comparisons
function [g1mean,g2mean,PVAL,COHEND,COHENLOWER,COHENUPPER] = statcompare(group1,group2,NAME1,NAME2,TESTTYPE)
%%% This function performs statistical testing to compare two groups (unpaired).

    if contains(TESTTYPE,"median")
        TEST = "mediandiff";
        YSTRING = "Median difference";
    else
        TEST = "cohen";
        YSTRING = "Cohen's d";
    end
    if isempty(NAME1)
        NAME1 = 'Group 1';
    end
    if isempty(NAME2)
        NAME2 = 'Group 2';
    end
    % Ensure column vectors
    if length(group1) == width(group1)
        group1 = group1.';
    end
    if length(group2) == width(group2)
        group2 = group2.';
    end
    % Combine into a single array
    if length(group1) > length(group2)
        stack = NaN(length(group1),2);
    else
        stack = NaN(length(group2),2);
    end
    stack(1:length(group1),1) = group1; stack(1:length(group2),2) = group2;
    % One-way ANOVA
    [PVAL,~,stats] = anova1(stack,[],'off');
    g1mean = stats.means(1);
    g2mean = stats.means(2);
    % Effect size
    effsize = meanEffectSize(group1,group2,Effect=TEST,...
        ConfidenceIntervalType="bootstrap", BootstrapOptions=statset...
        (UseParallel=true),NumBootstraps=3000);
    effsizearray = table2array(effsize);
    COHEND = effsizearray(1);
    COHENLOWER = effsizearray(2);
    COHENUPPER = effsizearray(3);
    ann(1) = {"p = "+string(PVAL)};
    ann(2) = {YSTRING+" = "+string(COHEND)+" ["+string(COHENLOWER)+...
        ", "+string(COHENUPPER)+"]"};

    gardnerAltmanPlot(group1,group2,Effect=TEST,...
    ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=3000);
    xticklabels({NAME1,NAME2,'Effect size'})
    yyaxis left; ylabel('Lifetime (ps)'); 
    yyaxis right; ylabel(YSTRING);
    annotation('textbox',[0.4 0.7 0.2 0.2],'String',ann,...
        'FitBoxToText','on','FontSize',12,'EdgeColor',[1 1 1]);
    title('')
    ax = gca; ax.FontSize = 12;
    return
end