from eewsimpy.util import inv2coord, haversine, tttinterp
from numpy import pi,argsort,nan,mgrid,zeros,nanmin,nanmax
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth

def leadtimes(inventory,
              catalog,
              target,
              declust=0.01, #degrees
              min_station_number=4,#minimum number of required station
              flat_latency={'SV.*.00.HN*':1.,
                            'GT*':5,
                            'SV*':2,
                            '*':3},
              mtodeg=2*6371000*pi/360,
              model='iasp91',
              debug=False
              ):
    
    """Compute Earthquake Early Warning (EEW) lead times.

    Each EEW lead time is defined for a given seismic event and EEW target. It is assumed as the delay between the arrival of the P-wave at the 4th closest station to the event and the arrival of the S-wave at the EEW target. However, in practice, EEW is required before the shaking at the EEW target exceeds the threshold for damage, which can be different from the S-wave arrival.

    .. code:: python

        from obspy.clients.fdsn.client import Client
        inv = Client('ETH').get_stations(level='channel')
        cat = Client('ETH').get_events(limit=1,orderby='magnitude',maxdepth=10000)
        from eewsimpy.event import leadtime
        leadtimes(inv, cat, target=[8.54690,47.37850])

    :param inventory: The instrument metadata inventory.
    :type inventory: :external:py:class:`obspy.core.inventory.inventory.Inventory` 
    :param catalog: The seismic event catalog.
    :type catalog: :external:py:class:`obspy.core.event.Catalog`
    :param target: The EEW target coordinates (longitude, latitude).
    :type target: :py:class:`list` of :py:class:`float`
    :param declust: Clustering threshold in degrees for station coordinates. Default is 0.01 degrees.
    :type declust: :py:class:`float`
    :param min_station_number: Minimum number of required stations. Default is 4.
    :type min_station_number: :py:class:`int`
    :param flat_latency: Dictionary specifying the flat latency values for different station patterns. Default values are provided.
    :type flat_latency: :py:class:`dict`
    :param mtodeg: Conversion factor from kilometers to degrees. Default is calculated based on Earth's radius.
    :type mtodeg: :py:class:`float`
    :param model: Time travel tables used for P-wave and S-wave travel times. Default is provided.
    :type model: :py:class:`str`
    :param debug: A boolean flag to enable/disable debug output. The default value is False.
    :type debug: :py:class:`bool`

    :return: EEW lead time (in seconds) for each event in the catalog.
    :rtype: :py:class:`numpy.ndarray`
    """

    if inventory is tuple and len(inventory)==3:
        statlats,statlons,statoff = inventory
    else:
        # Stations coordinates, one station per clusters
        statlats,statlons,statoff = inv2coord(inventory,
                                            declust=declust,
                                            flat_latency=flat_latency)    

    mod = TauPyModel(model=model)
    leadtime = zeros(len(catalog))
    for e,event in enumerate(catalog):

        origin = event.preferred_origin_id.get_referred_object()

        # Distance from event location to minimal number of stations 
        inputs=[origin.longitude,
                origin.latitude,
                statlons,
                statlats]
        dstations = haversine(*inputs)
        station = argsort(dstations)[min_station_number-1]
        inputs=[origin.longitude,
                origin.latitude,
                statlons[station],
                statlats[station]]
        dstation = gps2dist_azimuth(*inputs)[0]

        # Distance from event location to target 
        inputs=[origin.longitude,
                origin.latitude]
        dtarget = gps2dist_azimuth(*inputs,
                                   *target[:2])[0]
        
        
        # Delay to detect event considering its location
        opt={'source_depth_in_km':origin.depth/1000,
             'distance_in_degree':dstation/mtodeg,
             'phase_list':['ttp']} 
        ttp = min([a.time for a in mod.get_travel_times(**opt)])
        
        
        
        # EEW lead time between event detection and S-wave arrival at target
        opt['phase_list'] = ['tts']
        opt['distance_in_degree'] = dtarget/mtodeg
        tts = min([a.time for a in mod.get_travel_times(**opt)])

        if debug:
                print(event.short_str())
                print('ttS @ target: %.1f km in %s s & ttP @ %dth station: %.1f km in %s s'%(dtarget/1000,tts,min_station_number,dstation/1000,ttp))   

        leadtime[e] = tts - ttp

    return leadtime



def leadtimegrid(inventory,
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