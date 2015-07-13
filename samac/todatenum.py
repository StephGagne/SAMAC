################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####


def todatenum(d):
    # converts to ordinal number (days since jan 1 of year 1)
    # d is a datetime.datetime or a datetime.time object
    try: d1=d.toordinal()
    except: d1=0
    d2=(d.hour+(d.minute+d.second/60.)/60.)/24.
    dn=d1+d2
    return dn
    

