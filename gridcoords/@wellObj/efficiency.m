function o = efficiency(o,iComp)
    %WELLOBJ/EFFICIENCY -- computes and plots seasonal thermal efficiency
    % and ads the result to the  UserData field of the wellObj.
    %
    % well = well.efficiency(iComp)
    % iComp must be the species' index of the temperature in the simulation
    % and well must have the temperature on board in well(iw).Cout(iComp,:)
    % well = well.setCout(T,iComp) must, therefore, have already been
    % excecuted
    %
    % TO 120913
    
    aYear = 365.24;
    aWeek = 7;       % space around a year
    
    if isempty(o(1).Cout) || all(isnan(o(1).Cout(iComp,:)))
        error('%s: Cout not set, apply well=well.setCout(C,iComp) first',mfilename);
    end
    
    % Plot output concentration of wells versus time    
    for iw=1:length(o)
        I = o(iw).Q<0; % extraction episodes
        Eout =    I .* abs(o(iw).Q.*(o(iw).Cout(iComp,:)-o(iw).UserData.taverage));
        Ein  =  (~I).* abs(o(iw).Q.*(o(iw).C(   iComp,:)-o(iw).UserData.taverage));
        
        % get distance between times that are a year apart
        iEnd   = 1:length(o(iw).Dt);
        iStart = NaN(size(iEnd));
        for it=iEnd
            istrt = find(o(iw).t>=o(iw).t(it)-aYear-aWeek & o(iw).t<o(iw).t(it)-aYear+aWeek);
            if ~isempty(istrt),
                iStart(it)=istrt+1;
            end
        end
        if all(isnan(iStart))
            msgId = 'mfLab:seasonalATESwellEfficiency:simPeriodSmallerThanAYear';
            warning('on',msgId);
            beep;
            warning(msgId,...
                ['Time series for wells must be > 1 year to be able to compute its seasonal efficiency!\n',...
                'Remedy: make sure the total simulation time in the PER sheet is > 365 days, rather > 730 d']);
            warning('off',msgId);
        end        
        o(iw).UserData.efficiency = NaN(size(o(iw).Dt));
        for it=1:length(iStart)
            if ~isnan(iStart(it))
                o(iw).UserData.efficiency(it)=sum(Eout(iStart(it):iEnd(it)))./sum(Ein(iStart(it):iEnd(it)));
            end
        end
    end
end
