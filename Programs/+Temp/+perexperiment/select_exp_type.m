function[Dates, Fish, exptype, lum_th, intmax] = select_exp_type(char)

% ------- exp 60 -------
if char == 'exp60'
    disp('Chosen : EXP 60')
    Dates = {'17-11-14', '17-11-15','17-11-16', '17-11-29'};
    Fish = {'fish1', 'fish2', 'fish3', 'fish4', [];...
        'fish1', 'fish2', 'fish3', 'fish4', 'fish5';...
        'fish1', [], [], [], [];...
        'fish1', 'fish2', 'fish3', 'fish4', [] };
    
    exptype = 'utilisables_lum_exp';
    intmax = '60';
    lum_th = luminosity_exponentielle(str2double(intmax)/100);
    
% ------- exp 30 -------
elseif char == 'exp30'
    disp('Chosen : EXP 30')
    Dates = {'17-11-09', '17-11-16','17-11-21'};
    Fish = {'fish2', 'fish3', 'fish4', 'fish5', [];...
        'fish2', 'fish3', 'fish4', [], [];...
        'fish1', 'fish2', 'fish3', 'fish4', [] };
    
    exptype = 'utilisables_lum_exp';
    intmax = '30';
    lum_th = luminosity_exponentielle( str2double( intmax )/100);
    
% ------- sin 60 -------
elseif char == 'sin60'
    disp('Chosen : SIN 60')    
    Dates = {'17-11-22', '17-11-23','17-11-28'};
    Fish = {'fish1', 'fish2', 'fish3', 'fish4', 'fish5';...
        'fish1', 'fish2', 'fish3', 'fish4', [];...
        [], 'fish2', 'fish3', 'fish4', [] };
    
    exptype = 'utilisables_lum_sinus';
    intmax = '60';
    lum_th = luminosity_sinus( str2double( intmax )/100);
end

