from fnmatch import fnmatch
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
