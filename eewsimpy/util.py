from obspy.geodetics.base import gps2dist_azimuth
import numpy
from fnmatch import fnmatch
try:
    from shapely.geometry import LineString
    from shapely.ops import nearest_points
except Exception as e:
    print('shapely skipped')
    print(e)
try:
    import cartopy.io.shapereader as shpreader
except Exception as e:
    print('cartopy skipped')
    print(e)


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


from obspy.taup import TauPyModel
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

        from eewsimpy.util import eventcountrydistance
        events = [[ e.preferred_origin().longitude, 
                    e.preferred_origin().latitude, 
                    e.preferred_origin().depth/1000 ] for e in cat]
        countrydists = eventcountrydistance(country, events)

    :param countryname: The country name.
    :type countryname: :py:class:`string` 
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
