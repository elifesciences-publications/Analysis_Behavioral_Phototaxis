function[myfit] = AutocorrelationFit(data_to_fit_on, x, y, pturn, wturn, wfor, fig)
%% from data :

if size(x,1)<size(x,2)
    x=x';
end
if size(y,1)<size(y,2)
    y=y';
end

if sum(isnan(y))>0
    seqend = find(isnan(y), 1)-1;
    y = y(1 : seqend);
    x = x(1 : seqend);
end

switch data_to_fit_on
    case 1
        %% Fit : dX = f(dXn-1)
        %.........................................................................
        
        myfittype = fittype('autocorr_lag1_fit(x, p_flip, pturn, wturn, wfor)',...
            'problem', {'pturn', 'wturn', 'wfor'}, 'coefficients', {'p_flip'});
        
        myfit = fit(x,y,myfittype,'StartPoint', [0.2],'problem', {pturn, wturn, wfor});
        
        if fig
            figure
            plot(myfit, x, y)
        end
        
    case 2
        %% Fit : AC
        %.........................................................................
        
        myfittype = fittype('autocorr_lag10_fit(x, p_flip, pturn, wturn, wfor)',...
            'problem', {'pturn', 'wturn', 'wfor'}, 'coefficients', {'p_flip'});
        
        myfit = fit(x,y,myfittype,'StartPoint', [0.2],'problem', {pturn, wturn, wfor});
        
        if fig
            figure;
            plot(myfit, x, y)
        end
end