function[lum, dlum, x, dx, r, tb, ibi, t, f] = choose_experiment(e)

Temp.perexperiment.loadPooledExpsCell

if strcmp(e, 'exp60') % CHOOSE EXPERIMENT
    lum = e6_lum;
    dlum = e6_dlum;
    x = e6_x;
    dx = e6_dx;
    r = e6_r;
    t = e6_t;
    tb = e6_tb;
    ibi = diff(tb,1,2);
    f = e6_f;
elseif strcmp(e, 'exp30')
    lum = e3_lum;
    dlum = e3_dlum;
    x = e3_x;
    dx = e3_dx;
    r = e3_r;
    t = e3_t;
    tb = e3_tb;
    ibi = diff(tb,1,2);
    f = e3_f;
elseif strcmp(e, 'sin60')
    lum = s6_lum;
    dlum = s6_dlum;
    x = s6_x;
    dx = s6_dx;
    r = s6_r;
    t = s6_t;
    tb = s6_tb;
    ibi = diff(tb,1,2);
    f = s6_f;
end