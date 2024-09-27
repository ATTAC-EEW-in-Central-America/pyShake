# -*- coding: utf-8 -*-
"""
Misc utility module. 
"""

from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event import Catalog
from obspy import read_events, UTCDateTime

from os.path import exists

import subprocess

import numpy
from fnmatch import fnmatch
try:
    from shapely.geometry import LineString, Polygon
    from shapely.ops import nearest_points
except Exception as e:
    print('shapely skipped')
    print(e)
try:
    import cartopy.io.shapereader as shpreader
except Exception as e:
    print('cartopy skipped')
    print(e)

from pyshake import gm



def inv2coord(inventory,
              declust=0.01,
              flat_latency={'*':0}):
    """
    Returns the list of stations latitudes and longitudes, 
    omiting old stations being too close to another one
    """
    statlats=[s.latitude for n in inventory for s in n ]
    statlons=[s.longitude  for n in inventory for s in n ]
    statoff=[]
    statid=[]
    for n in inventory:
        for s in n :
            for c in s :
                mem=0
                mseedid = '%s.%s.%s.%s'%(n.code,s.code,c.location_code,c.code)
                for k in flat_latency:
                    if (len(k)>mem and 
                        fnmatch(mseedid,k)):
                        mem=len(k)
                        delay=flat_latency[k]
                        memid=mseedid
            statoff+=[delay]
            statid+=[mseedid]
            
    tooclose=[]
    statindexes=range(len(statlats))
    for i in statindexes:
        for j in range(i+1,len(statlats)):
            if (abs(statlats[i]-statlats[j])<declust and 
                abs(statlons[i]-statlons[j])<declust):
                tooclose+=[i]
                if statoff[i]>statoff[j]:
                    statoff[i]=statoff[j]
    statlats = [statlats[i] for i in statindexes if i not in tooclose]
    statlons = [statlons[i] for i in statindexes if i not in tooclose]
    statoff = [statoff[i] for i in statindexes if i not in tooclose]

    return statlats,statlons,statoff


from numpy import deg2rad,sin,cos,sin,arcsin,sqrt,abs
def haversine(lon1, lat1, lon2, lat2,
              r = 6371 # Radius of earth in kilometers. Use 3956 for miles
              ):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(deg2rad, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * arcsin(sqrt(a))
    
    return abs(c * r)

taup=True
try:
    from obspy.taup import TauPyModel
except Exception as e:
    print(e)
    taup=False
from scipy.interpolate import interp1d
from numpy import linspace
def tttinterp(model='iasp91',
              depth=5,
              ttt={'p':None,'s':None},
              dmax=3):
    """
    Return the dict of interpolants for p and s-phase 
    travel times modeling at any distance such as:
    `time = ttt[phases](dist)`
    """
    if not taup:
        print('TauP cannot be imported, this does not work')
        return ttt
    mod = TauPyModel(model=model)
    distances=linspace(.001,dmax,64)
    
    for phases in ttt.keys():
        arrivals=[]
        for d in distances:
            opt={'source_depth_in_km':depth,
                 'distance_in_degree':d,
                 'phase_list':['tt'+phases],
                 'receiver_depth_in_km':0.0}
            arrivals+=[min([a.time for a in mod.get_travel_times(**opt)])]
        ttt[phases] = interp1d(distances,arrivals)
        #time = ttt[phases](dist)
    return ttt

                
from geopip import search
def isincountries(la, lo, countrycodes,
                  dthresh=1.4):
    """
    The function `isincountries` checks if a given latitude and longitude pair 
    is within a specified distance threshold of any country in a list of ISO2
    country codes. The function uses the `geopip` library to determine the country 
    code for the given coordinates
    
    :param la: Input latitude.
    :type la: :py:class:`float`
    :param lo:  Input longitude.
    :type lo: :py:class:`float`
    :param countrycodes: A comma separated list of ISO2 country codes.
    :type countrycodes: :py:class:`str`
    :param dthresh: The distance threshold within which the function will search 
        for countries.
    :type dthresh: :py:class:`float`
    
    :return: The function will return `True` if the coordinates are within any of the 
        specified countries, and `False` otherwise.
    :rtype: :py:class:`bool`
    """

    if countrycodes is None:
        return True
    inside=False
    lalo=[la,lo]
    for lad in [la+dthresh,la-dthresh,la]:
        for lod in [lo+dthresh,lo-dthresh,lo]:
            gcode=search(lng=lod, lat=lad)
            if gcode is None or gcode['ISO2'] is None:
                continue
            if gcode['ISO2'].lower() in countrycodes:
                return True
    if gcode is None or gcode['ISO2'] is None:#gcode.country_code is None:
        return
    

def eventdistance(la, lo, de, 
                  refla, reflo, refde,
                  v=1,
                  meters=800):
    """
    The function `eventdistance` calculates the distance between two geographical 
    coordinates, taking into account the depth of the locations if provided.
    
    :param la: Latitude of the first location.
    :param la: :py:class:`float`
    :param lo: Longitude of the first location.
    :param lo: :py:class:`float`
    :param de: Depth of the first location (set to None if unknown).
    :param de: :py:class:`float`
    :param refla: Latitude of the second location.
    :param refla: :py:class:`float`
    :param reflo: Longitude of the second location.
    :param reflo: :py:class:`float`
    :param refde: Depth of the second location (set to None if unknown).
    :param refde: :py:class:`float`
    :param v: Divider of the horizontal distance output in case depth is unknown.
    :param v: :py:class:`float`
    :param meters: If the absolute value of the depth (`de` or `refde`) is 
        greater than this threshold, the function assumes that the depth is 
        in kilometers instead of meters. 
    :param meters: :py:class:`float`
    
    :return: The distance in kilometers between the two geographical coordinates.
    :rtype: :py:class:`float`
    """
    laloref=[la,lo]
    lalo = [refla, reflo]
    #d = 111.2*locations2degrees(lalo[0],lalo[1],laloref[0],laloref[1])
    d = gps2dist_azimuth(lalo[0],lalo[1],laloref[0],laloref[1])[0]/1000

    if refde is None or de is None:
        return d/v

    if abs(de)>meters:
        print('Assuming depth %s in m'%de)
        de = de/1000
    if abs(refde)>meters:
        print('Assuming depth %s in m'%refde)
        refde = refde/1000

    return  ((d**2 + (de-refde)**2)**.5)/v


def eventcountrydistance(countryname, 
                         events, 
                         res=1, 
                         defaultdepth=5,
                         defaultepd=2,
                         quickndirty=False):
    
    """Compute the distance from given hypocenter(s) to given country.

    In case the hypocenter is within the country, the output distance is taken from hypocenter depth (default to 5 km).

    .. code:: python

        from obspy.clients.fdsn.client import Client
        cat = Client('ETH').get_events(limit=1,orderby='magnitude',maxdepth=10000)

        from pyshake.util import eventcountrydistance
        events = [[ e.preferred_origin().longitude, 
                    e.preferred_origin().latitude, 
                    e.preferred_origin().depth/1000 ] for e in cat]
        countrydists = eventcountrydistance(country, events)

    :param countryname: The country name.
    :type countryname: :py:class:`str` 
    :param events: The event(s) hypocenter coordinates.
    :type events: :py:class:`list` of :py:class:`list` e.g., [[lo1,la1],[lo2,la2,z2]].

    :return: Country distance for each event in the input list.
    :rtype: :py:class:`list`
    """
    
    single=False
    if not isinstance(events[0],(list,tuple)) :
        single=True
        events=[events]

    if False in [len(e)>2 for e in events]:
        print('ERROR: Input(s) are not at least 2d:')
        print(events)

    for e,event in enumerate(events):
        if len(event)==5:
            elon = event[0] + event[3]/111.2 * numpy.sin(numpy.deg2rad(event[4]))
            elat = event[1] + event[3]/111.2 * numpy.cos(numpy.deg2rad(event[4]))

            elon = [event[0]]+[l for l in numpy.linspace(event[0],elon,int(event[3]/res))]+[elon]
            elat = [event[1]]+[l for l in numpy.linspace(event[1],elat,int(event[3]/res))]+[elat]

            events[e] = LineString([(elon[p],elat[p]) for p in range(len(elon))])
            events[e].depth = event[2]
        elif len(event)==3:
            events[e] = LineString([(event[0],event[1]), (event[0],event[1])])
            events[e].depth = event[2]
        elif len(event)==2:
            events[e] = LineString([(event[0],event[1]), (event[0],event[1])])
            events[e].depth = defaultdepth
        else:
            print('ERROR: Input are not at least 2d:')
            print(event)
        if events[e].depth is None:
            print('WARNING: Event (%s) depth is None, using %s default depth)'%(','.join(['%s'%f for f in event]), defaultdepth))
            events[e].depth = defaultdepth
    
    distances=[None for e in events]
    url = shpreader.natural_earth(resolution='110m',
                                        category='cultural',
                                        name='admin_0_countries')
    reader = shpreader.Reader(url)
    countries = reader.records()
    country = [country for country in countries        if countryname.lower() in country.attributes['NAME'].lower() ]
    if not len(country):
        print('ERROR: Could not find %s in [%s]'%(countryname,','.join([country.attributes['NAME'] for country in countries])))
    else:
        if len(country)>1:
            print('WARNING: Considering the first of [%s] as matching [%s]'%(','.join([c.attributes['NAME'] for c in country]),countryname))
        country=country[0]
        
        for e,event in enumerate(events):
            if country.geometry.intersects(event):
                #print('INFO: Event intersects %s using depth %d km as country distance'%(countryname, event.depth))
                distances[e] = ( event.depth**2 + defaultepd**2 )**.5
            else:
                if quickndirty:
                    print('WARNING: Event outside %s using shapely distance as country distance (can be VERY incorrect)'%(countryname))
                    distances[e] = ((country.geometry.distance(event)*111.2)**2 + event.depth**2)**.5 
                else:
                    impact_trajectory = nearest_points(event,country.geometry)
                    impact_trajectory = [c for p in impact_trajectory for c in p.coords[:][0][::-1]]
                    distances[e] = ((gps2dist_azimuth(*impact_trajectory)[0]/1000)**2 + event.depth**2)**.5
    if single:
        return distances[0]
    return distances


    
def polygonthreshold(polygon,
                     countryname='Switzerland',
                     minintensity=4,
                     mineventdepth=1,
                     minepicentraldistance=1,
                     groundmotionmodel=gm.gmm,
                     quickndirty=False,
                     **kwargs):
    """
    The function `polygonthreshold` calculates the minimum magnitude of an earthquake 
    required to reach a given shaking intensity level anywhere inside a given country.

    .. code:: python

        from pyshake.util import polygonthreshold

        # a polygon over Zurich
        polygon=[(8,47), 
                (9,47), 
                (9,47.5), 
                (8,47.5)]

        # the minimum magnitude to reach MMI 3 in Switzerland 
        # from any event in the polygon over Zurich
        # `s = 0` : No amplification
        polygonthreshold(polygon, 
                        countryname = 'Switzerland', 
                        minintensity = 3, 
                        s = 0) 
        # 0.97

        # the minimum magnitude to reach MMI 3 in Italy 
        # from any event in the polygon over Zurich
        # `s = 1` : An amplification of 1 intensity unit 
        polygonthreshold(polygon, 
                        countryname = 'Italy', 
                        minintensity = 3, 
                        s = 1) 
        # 3.99

    :param polygon: The polygon parameter is a set of coordinates that define the shape 
        of the polygon. It is used to determine if the event intersects with the specified 
        country
    :type polygon: :py:class:`list` of :py:class:`list` e.g., [[lo1,la1],[lo2,la2]].
    :param countryname: The name of the country where the intensity is evaluated. The 
        default value is 'Switzerland', defaults to Switzerland (optional)
    :type countryname: :py:class:`str` 
    :param minintensity: The minimum intensity threshold for the polygon. Any earthquake 
        with an intensity below this threshold will not be considered, defaults to 4 
        (optional)
    :type minintensity: :py:class:`float`    
    :param mineventdepth: The mineventdepth parameter represents the depth of the event 
        in kilometers. It is used to calculate the distance between the event and the 
        polygon, defaults to 1 (optional)
    :type mineventdepth: :py:class:`float`
    :param minepicentraldistance: The parameter "minepicentraldistance" represents the 
        minimum epicentral distance in kilometers. It is used in the calculation of the 
        distance between the event and the country, defaults to 1 (optional)
    :type minepicentraldistance: :py:class:`float`
    :param groundmotionmodel: The groundmotionmodel parameter is a function that 
        calculates the intensity of ground motion given the distance and magnitude of an 
        earthquake. In this code, it is set to gmm.gm, which is a reference to a default 
        ground motion model function by Allen et al. (2012)
    :type groundmotionmodel: :class:`pyshake.gm.gmm`
    :param quickndirty: The parameter "quickndirty" is a boolean flag that determines 
        whether a quick and dirty calculation should be used for the distance calculation 
        if the event is outside the specified country. If set to True, the shapely 
        distance between the event and the country will be used, which can be very 
        incorrect, defaults to False (optional)
    :type quickndirty: :py:class:`bool`

    :return: the minimum magnitude of an earthquake that would produce a ground motion 
        intensity greater than or equal to the specified minimum intensity, given the 
        parameters and inputs provided.
    :rtype: :py:class:`float`
    """

    polygon = Polygon(polygon)
    #polygon.depth = mineventdepth
    
    url = shpreader.natural_earth(resolution='110m',
                                        category='cultural',
                                        name='admin_0_countries')
    reader = shpreader.Reader(url)
    countries = reader.records()
    country = [country for country in countries        if countryname.lower() in country.attributes['NAME'].lower() ]
    if not len(country):
        print('ERROR: Could not find %s in [%s]'%(countryname,','.join([country.attributes['NAME'] for country in countries])))
    else:
        if len(country)>1:
            print('WARNING: Considering the first of [%s] as matching [%s]'%(','.join([c.attributes['NAME'] for c in country]),countryname))
        country=country[0]
        
        if country.geometry.intersects(polygon):
            print('INFO: Event intersects %s using depth %.1f km and  minimal epicentral distance of %.1f km as distance'%(countryname, mineventdepth, minepicentraldistance))
            distance = ( mineventdepth**2 + minepicentraldistance**2 )**.5
        else:
            if quickndirty:
                print('WARNING: Event outside %s using shapely distance as country distance (can be VERY incorrect)'%(countryname))
                distance = ((country.geometry.distance(polygon)*111.2)**2 + mineventdepth**2)**.5 
            else:
                impact_trajectory = nearest_points(polygon,country.geometry)
                impact_trajectory = [c for p in impact_trajectory for c in p.coords[:][0][::-1]]
                distance = ((gps2dist_azimuth(*impact_trajectory)[0]/1000)**2 + mineventdepth**2)**.5
    intensity = 99
    for minmag in numpy.arange(10,-10,-0.01):
        if intensity <= minintensity:
            return minmag
        intensity = groundmotionmodel.get_intensity(distance,minmag, **kwargs)
    return None

def ssh_event(sshuserhost='user@host',
              dbplugin='dbmysql',
              db='localhost/seiscomp3',
              opt='-fma',
              seiscomp='',
              eventid=None,
              minlatitude=-90,
              maxlatitude=90,
              minlongitude=-180,
              maxlongitude=180,
              minmagnitude=0,
              maxmagnitude=10,
              mintime=UTCDateTime()-60*60,
              maxtime=UTCDateTime()):

    ssh = subprocess.Popen(["ssh",
                            "-i ~/.ssh/id_rsa",
                            sshuserhost],
                            stdin =subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True,
                            bufsize=0)
    
    if eventid is None:
        queryparam = (minlatitude,
                    maxlatitude,
                    minlongitude,
                    maxlongitude,
                    minmagnitude,
                    maxmagnitude,
                    mintime.isoformat().replace('T',' ')[:19],
                    maxtime.isoformat().replace('T',' ')[:19])
        query="SELECT PEvent.publicID, Origin.time_value AS OT, Origin.latitude_value,Origin.longitude_value, Origin.depth_value, Magnitude.magnitude_value, Magnitude.type FROM Origin,PublicObject as POrigin, Event, PublicObject AS PEvent, Magnitude, PublicObject as PMagnitude WHERE Event._oid = PEvent._oid AND Origin._oid = POrigin._oid AND Magnitude._oid = PMagnitude._oid AND PMagnitude.publicID=Event.preferredMagnitudeID AND POrigin.publicID = Event.preferredOriginID AND Origin.latitude_value >= %s AND Origin.latitude_value <= %s AND Origin.longitude_value >= %s AND Origin.longitude_value <= %s AND Magnitude.magnitude_value >= %s AND Magnitude.magnitude_value <= %s AND Origin.time_value >= '%s' AND Origin.time_value <= '%s';"%queryparam
        scquery = '%s scquery --debug --plugins %s -d %s -Q \"%s\" |grep -v "^$"| while read ID TRASH;do %s scxmldump %s --debug --plugins %s -d %s   -E $ID ;done' % (seiscomp, dbplugin, db, query, seiscomp, opt, dbplugin, db)
    
        ssh.stdin.write(scquery)
    else :
        #%s scquery --debug --plugins %s -d %s -Q \"%s\" |grep -v "^$"| while read ID TRASH;do % seiscomp, dbplugin, db, query, 
        scdump = '%s scxmldump %s --debug --plugins %s -d %s -E %s ' % (seiscomp, opt, dbplugin, db, eventid)

        ssh.stdin.write(scdump)

    ssh.stdin.close()

    # Fetch output
    string = ''.join([line.replace('>\n','>') for line in ssh.stdout if ('<' in line and '>' in line)])
    if len(string):
        string=string.replace('0.11" version="0.11','0.10" version="0.10')
        string=string.replace('0.12" version="0.12','0.10" version="0.10')
        f=open('/tmp/tmp.xml', 'w+') 
        f.write(string)
        f.close()
        try:
            return read_events('/tmp/tmp.xml',format='SC3ML')
        except:
            for n in range(9999):
                path_to_file='invalid.%s.xml'%n
                if not exists(path_to_file):
                    break
            f=open(path_to_file, 'w+') 
            f.write(string)
            f.close()
            print('Cannot read',path_to_file)
            return Catalog()
    else:
        #print('Failed:')
        #print(scquery)
        #print([line for line in ssh.stderr])
        string = ''.join([line for line in ssh.stderr])
        f=open('.stderr', 'w+') 
        f.write(string)
        f.close()
        print('Cannot get event from %s .stderr'%sshuserhost)
        return Catalog()

