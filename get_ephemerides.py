import urllib
import requests

code = {'Mercury':'199', 'Venus':'299', 'Earth':'399', 'Mars':'499',
        'Jupiter':'599', 'Io':'501', 'Europa':'502', 'Ganymede':'503',
        'Saturn':'699', 'Uranus':'799', 'Neptune':'899','Callisto':'504',
        'EUROPA': '502'}
# add the target name to this list to use the fits header 'TARGNAME'

def get_ephemerides(target,Tstart,Tend) :

    dt      = 22*3600. # time interval [s]
    http = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1"
    make_ephem = "&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'"
    command    = "&COMMAND=" + code[ target ]
    center     = "&CENTER='@-170'"                            # HST = '@-48', ALMA = '-7', JWST = '-170'
    t_start = "&START_TIME='" + Tstart + "'"
    t_stop = "&STOP_TIME='" + Tend + "'"
    t_step     = "&STEP_SIZE='1m'"
    quantities = "&QUANTITIES='13,14,17,15,19,20'" #13 = angular diameter in arcseconds, 14 = sub. obs. lon and lat, 15 = sub. sun lon and lat, 17 = np. angle and angular distance from center of disk, 19 is heliocentric range and range rate in AU and km/s, 20 is observer range and range rate.
    csv        = "&CSV_FORMAT='YES'"

    url = http+make_ephem+command+center+t_start+t_stop+t_step+quantities+csv

    ephem = requests.get(url)
    ephem = ephem.text.splitlines(True)

    data = []
    for i in range(len(ephem)) :
        if ephem[i].startswith("$$SOE") :
            j = i+1
            while not ephem[j].startswith("$$EOE") :
                data.append(ephem[j].split(','))
                j+=1
    return data
