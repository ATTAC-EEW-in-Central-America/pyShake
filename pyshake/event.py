from pyshake.util import inv2coord, haversine, tttinterp
from numpy import pi,argsort,nan,mgrid,zeros,nanmin,nanmax
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth

def leadtimes(inventory,
              catalog,
              target,
              declust=0.01, 
              min_station_number=4,
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
        from pyshake.event import leadtime
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
                                              declust=declust)    

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
    """
    Calculates the EEW delays for a grid of event locations based on station coordinates and other parameters.

    :param inventory: The instrument metadata inventory, which contains information about the seismic stations in the network
    :type inventory: :external:py:class:`obspy.core.inventory.inventory.Inventory`

    :param declust: The declust parameter is the clustering distance in degrees. It determines the maximum distance between two stations for them to be considered part of the same cluster. Stations within the clustering distance will be grouped together and only one station from each cluster will be used in the calculations
    
    :param min_station_number: The `min_station_number` parameter is the minimum number of required stations for the calculation. It specifies the minimum number of stations that need to detect an earthquake in order for the calculation to be performed. If the number of stations that detect the earthquake is less than `min_station_number`, the calculation will, defaults to 4 (optional)
    
    :param target: The target parameter is the coordinates of the target location where you want to calculate the EEW delays. In this case, the target location is set to Guatemala City with latitude -90.535278 and longitude 14.613333
    
    :param depth: The depth parameter represents the depth of the earthquake events in kilometers, defaults to 5 (optional)
    
    :param dlat: The parameter `dlat` represents the latitude increment used to create a grid of latitude values. It determines the spacing between adjacent latitude values in the grid 
    
    :param dlon: The parameter `dlon` represents the longitude increment used to create a grid of longitude values. It determines the spacing between adjacent longitude values in the grid
    
    :param dmax: The parameter `dmax` represents the maximum distance in degrees from the network that will be considered for calculating the EEW delays, defaults to 3 (optional) 
    
    :param flat_latency: The `flat_latency` parameter is a dictionary that specifies the flat latency values for different station names or patterns. The keys in the dictionary are station names or patterns, and the values are the corresponding flat latency values. The flat latency is the delay in seconds that is added to the EEW delay for 
    
    :param resampling: The `resampling` parameter determines the resampling factor for the EEW delays. It controls the spacing between the grid points in the output EEW delays. A value of 1 means no resampling, while a value greater than 1 will result in a coarser grid with larger spacing between, defaults to 1 (optional)
    
    :param kmtodeg: The parameter "kmtodeg" is a conversion factor used to convert distances in kilometers to degrees. It is calculated as 2 times the product of the Earth's radius (6371 km) and pi divided by 360. This conversion factor is used in the calculations of distances between stations and source grid nodes
    
    :param Vp: Vp is the P-wave velocity, which is the speed at which P-waves travel through the Earth's crust. It is used to calculate the travel time of seismic waves from the event location to the target location
    
    :param ttt: The parameter `ttt` is an object of the `tttinterp` class, which is used for interpolating travel times. 
    :rtype: :py:class:`function`
    
    :return: three values: longitudes, latitudes, and EEW delays.
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