from fnmatch import fnmatch
def inv2coord(inventory,
              declust=0.01,
              flat_latency={'*':3}):
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


from numpy import deg2rad,sin,cos,sin,arcsin,sqrt
def haversine(lon1, lat1, lon2, lat2,
              r = 6371 # Radius of earth in meters. Use 3956 for miles
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
    
    return c * r


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


from numpy import pi,sort,nan,mgrid,zeros,nanmin,nanmax
def eewtargetgrid(inventory,
                  declust=0.01, #degrees
                  min_station_number=4,#minimum number of required station
                  target=[-90.535278,14.613333],#Guatemala City
                  depth=5,#of event in km 
                  dlat=.1,#degrees
                  dlon=.1,#degrees
                  dmax=3,#degrees
                  flat_latency={'SV.*.00.HN*':1.,
                                'GT-rsn2.*':1.,
                                'GT-altac.*':1.,
                                'GT*':5,
                                'SV*':2,
                                '*':3},
                  resampling=1,
                  kmtodeg=2*6371*pi/360,
                  Vp=5.5,
                  ttt = tttinterp(),
                  ):
    
    """...

    :param inventory: The instrument metadata inventory.
    :type inventory: :py:class:`obspy.core.inventory.inventory`
    :returns: Longitudes, latitudes, and EEW delays
    :rtype: :py:class:`tuple`
    """
    
    # Stations coordinates, one station per clusters
    statlats,statlons,statoff = inv2coord(inventory,
                                          declust=declust,
                                          flat_latency=flat_latency)    
    print(len(statlons),'stations')
    
    # Grids within and around network
    i=slice(nanmin(statlats)-(dmax),
            nanmax(statlats)+(dmax),
            dlat)
    j=slice(nanmin(statlons)-(dmax),
            nanmax(statlons)+(dmax),
            dlon)
    latitudes, longitudes = mgrid[i,j]

    opt=[latitudes.shape[0],latitudes.shape[1]]
    EEWgrid = zeros(opt)
    
    # Iterating over all event locations
    for i in range(latitudes.shape[0]):
        for j in range(latitudes.shape[1]):

            # Distance to minimal station 
            # number from current event location
            inputs=[longitudes[i,j],
                    latitudes[i,j],
                    statlons,statlats]
            dstations = haversine(*inputs)
            dstations = [statoff[k]*Vp+d for k,d in enumerate(dstations)]
            dstation = sort(dstations)[min_station_number-1]
            
            # Distance to target from current event location
            inputs=[longitudes[i,j],latitudes[i,j]]
            dtarget= haversine(*inputs,
                               *target)   

            if (dstation>dmax*kmtodeg or 
                dtarget>dmax*kmtodeg):
                EEWgrid[i,j]=nan
                continue
            
            # Delay to detect event at this location
            EEWdelay = ttt['p'](dstation/kmtodeg)
            
            # EEW head time between event detection and S-wave arrival at target
            EEWgrid[i,j] = ttt['s'](dtarget/kmtodeg)-EEWdelay

    return longitudes,latitudes,EEWgrid

import time
def eewleadtime(inventory,
                catalog,
                target=[-90.535278,14.613333],#Guatemala City
                declust=0.01, #degrees
                min_station_number=4,#minimum number of required station
                depth=5,#of event in km 
                flat_latency={'SV.*.00.HN*':1.,
                              'GT*':5,
                              'SV*':2,
                              '*':3},
                resampling=1,
                kmtodeg=2*6371*pi/360,
                Vp=5.5,
                ttt = tttinterp()):
    
    """Compute EEW lead time

    The EEW lead time is defined for a given seismic event and EEW target. It is assumed as the delay between the arrival of the P wave at the 4th event-closest station and the arrival at S wave at the EEW target. However, in practice, EEW is required before the shaking at EEW target exceeds the threshold for damage which can be different than S-wave arrival.

    :example:

    >>> from obspy.clients.fdsn.client import Client
    >>> inv = Client('ETH').get_stations(level='channel')
    >>> cat = Client('ETH').get_events(limit=1,orderby='magnitude')
    >>> from eewsimpy.inv2tt import eewleadtime
    >>> eewleadtime(inv, cat, target=[8.54690,47.37850])
    array([17.04239037])

    :param inventory: The instrument metadata inventory.
    :type inventory: :py:class:`obspy.core.inventory.inventory`
    :param catalog: The seismic event catalog.
    :type catalog: :py:class:`obspy.core.catalog.catalog`
    :param target: The EEW target coordinates (e.i. longitude, latitude).
    :type target: :py:class:`list` of :py:class:`float`
    :returns: EEW lead time (in seconds) for each event in catalog.
    :rtype: :py:class:`list`
    """
    
    tic = time.perf_counter()
    if inventory is tuple and len(inventory)==3:
        statlats,statlons,statoff = inventory
    else:
        # Stations coordinates, one station per clusters
        statlats,statlons,statoff = inv2coord(inventory,
                                            declust=declust,
                                            flat_latency=flat_latency)    
    print(len(statlons),'stations')
    toc = time.perf_counter()
    print(f" {toc - tic:0.4f} seconds")

    leadtime = zeros(len(catalog))
    for e,event in enumerate(catalog):

        tic = time.perf_counter()
        print(f" {tic - toc:0.4f} seconds")
        # Distance from event location to minimal number of stations 
        inputs=[event.preferred_origin_id.get_referred_object().longitude,
                event.preferred_origin_id.get_referred_object().latitude,
                statlons,
                statlats]
        dstations = haversine(*inputs)
        dstations = [statoff[k]*Vp+d for k,d in enumerate(dstations)]
        dstation = sort(dstations)[min_station_number-1]
        
        toc = time.perf_counter()
        print(f" {toc - tic:0.4f} seconds")

        # Distance from event location to target 
        inputs=[event.preferred_origin_id.get_referred_object().longitude,
                event.preferred_origin_id.get_referred_object().latitude]
        dtarget= haversine(*inputs,
                           *target)   
        

        tic = time.perf_counter()
        print(f" {tic - toc:0.4f} seconds")

        # Delay to detect event considering its location
        EEWdelay = ttt['p'](dstation/kmtodeg)
        
        # EEW lead time between event detection and S-wave arrival at target
        leadtime[e] = ttt['s'](dtarget/kmtodeg)-EEWdelay

        toc = time.perf_counter()
        print(f" {toc - tic:0.4f} seconds")

    return leadtime

