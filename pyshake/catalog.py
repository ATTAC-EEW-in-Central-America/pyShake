# -*- coding: utf-8 -*-
"""
Catalog module. 
"""

import io, glob

from numpy import argmin, asarray, sqrt, argsort, nan, sort
from pandas import read_csv, DataFrame, concat
import cartopy.crs as ccrs

from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics.base import gps2dist_azimuth

from pyshake.util import isincountries, eventdistance, eventcountrydistance
from pyshake.plotting import plot, bmap, size2mag, mag2size
from pyshake.gm import gmm

import imp
import pyshake.plotting
imp.reload(pyshake.plotting)
plot = pyshake.plotting.plot

from matplotlib.pyplot import figure
from matplotlib import colors
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm 
from matplotlib.colors import NoNorm
from matplotlib.legend_handler import HandlerPathCollection
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#print(Tp)
#    0         1      2      3                    4                                5      6      7       8     9             10                          11
#[-91.3759, 13.988, 78.238, 5.4, UTCDateTime(2023, 12, 14, 17, 55, 31, 836000), -9.836, -91.58, 13.99, 62.51, 5.15, UTCDateTime(2023, 12, 14, 17, 55), 28.93],
#   lon       lat      dep   ma                 time                             bad dt   lon    lat     dep   ma         time                       ref dt
#    rf        ref      ref  ref                 ref                               ?      EEW    eew    eew   eew         eew                        eew                 
#print(Fp)
#print(Fn)

def sortarg(x):
    s = sorted(x)
    return [s.index(d) for d in x]

def get_eventindex(url,method='pandas',**opt):
    """
    Rapid event index list from FDSN web service. Uses the event 
    FDSN web-service, returns the list of events.

    .. code:: python

        from pyshake.catalog import get_eventindex
        get_eventindex('http://eida.ethz.ch', limit=100)

    :param url: The fdsnws base url (e.g. 'http://eida.ethz.ch').
    :type string: :py:class:`string`.
    :param method: 'pandas' or 'obspy' defines whether the event list 
        should returned as a :external:py:mod:`pandas.DataFrame` or a
        :external:py:mod:`obspy.core.event.Catalog`.
    :type string: :py:class:`string`.

    .. Note::

        Any additional keyword arguments will be passed to the web-service as
        additional arguments. Passing any non-default parameters that the 
        web-service does not support will raise an error.

    :return: The list of events.
    :rtype: :external:py:mod:`pandas.DataFrame` or :external:py:mod:`obspy.core.event.Catalog`
    """
    if method == 'pandas':

        url += '/fdsnws/event/1/query?format=text&'
        url += '&'.join(['%s=%s'%(k,opt[k])for k in opt])
        print(url)

        return read_csv(url, delimiter="|")
    
    elif method == 'obspy':

        return Client("IRIS").get_events(format='text',**opt)
    
    
def eew_report2dict(files='.seiscomp/EEW_reports/*txt'):
    """
    Reads earthquake early warning (EEW) report files from a specified 
    file-path pattern, and returns a dictionary of the extracted information. 
    It supports the following formats:

    .. code:: shell

        Tdiff |Type|Mag.|Lat.  |Lon.   |Depth |origin time (UTC)      |Lik.|Or.|Ma.|Str.|Len. |Author   |Creation t.            |Tdiff(current o.)
        ------------------------------------------------------------------------------------------------------------------------------------------
         80.60| MVS|3.42| 15.00| -92.89|  5.00|2019-09-24T20:21:01.06Z|0.99|  8|  4|    |     |scvsmag@s|2019-09-24T20:22:21.66Z| 80.60

        Mag.|Lat.  |Lon.   |tdiff |Depth |creation time (UTC)      |origin time (UTC)        |likeh.|#st.(org.) |#st.(mag.)
        ------------------------------------------------------------------------------------------------------------------
        5.15| 10.31| -86.27| 47.14| 10.00|2017-03-17T10:32:56.1282Z|2017-03-17T10:32:08.9848Z|  0.30|         19|         7
        
    :param files: A string that represents the file path pattern for the EEW report 
        files to be processed.
    :type files: :py:class:`str`
    
    :return: A dictionary where the keys are event identifiers (from file paths of EEW report 
        files from the specified file path pattern), and the values are dictionaries 
        representing the contents of the EEW reports.
    :rtype: :py:class:`dict`
    """
    
    reports={}
    for reportfile in glob.glob(files):
        with open(reportfile) as lines:
            lines = [[elt.replace(' ','')[:19] for elt in line.split('|')] for line in lines.read().splitlines() if '--' not in line and ' |#St.   | ' not in line]
            if not len(lines):
                print(reportfile, 'empty')
                continue
                reports[reportfile] = {}

            report = {k:[l[i] for l in lines[2:]] for i,k in enumerate(lines[0])}
            for k in report:
                if 'time' in k or 'reation' in k  or 'origin' in k:
                    report[k] = [UTCDateTime(e) for e in report[k]]  
                elif k in ['Type' , 'Author'] :
                    continue
                else:
                    try :
                        report[k] = [e.isspace() or not len(e) or float(e) for e in report[k]]  
                    except:
                        print('invalid!!!',reportfile)
                        print(lines[:3],k)

            for k in report:
                for i,e in enumerate(report[k]):
                    if isinstance(e,bool) and e==True:
                        report[k][i]=None

            if 'creation time (UTC)' in report:#'creationtime(UTC)', 'origintime(UTC)'
                report['Creation t.']=report['creation time (UTC)']
            if 'creationtime(UTC)' in report:#'creationtime(UTC)', 'origintime(UTC)'
                report['Creation t.']=report['creationtime(UTC)']
            if 'Creationt.' in report:
                report['Creation t.']=report['Creationt.']
            if 'origintime(UTC)' in report:#'creationtime(UTC)', 'origintime(UTC)'
                report['origin time (UTC)']=report['origintime(UTC)']
            if 'Type' not in report:
                report['Type']=['MVS' for e in report['origintime(UTC)']]
            if '#st.(org.)' in report:
                try:
                    report['Tdiff']=report['tdiff']
                    report['Lik.']=report['likeh.']
                    report['Or.']=report['#st.(org.)']
                    report['Ma.']=report['#st.(mag.)']
                    report['Str.']=[None for e in report['#st.(mag.)']]
                    report['Len.']=[None for e in report['#st.(mag.)']]
                except:
                    print(report.keys())
                    print('aliases undone')

            reports[reportfile] = report

    return reports

def reportfilter(reports,
                 countrycodes='ni',
                 dthresh=1.4
                 ):
    """
    Filters out reports based on specified country codes and distance
    threshold, removing duplicate locations and empty reports.
    
    :param reports: A dictionary of EEW reports as provided by 
        :func:`pyshake.catalog.eew_report2dict`.
    :type reports: :py:class:`dict`
    :param countrycodes: A string that specifies the country code or 
        codes to filter the reports by. By default, it is set to 'ni', 
        which represents the country code for Nicaragua. 
    :type countrycodes: :py:class:`str`, optional
    :param dthresh: The distance threshold (in degrees) used to determine 
        if a given (latitude and longitude) point is close enough to the 
        specified country to be kept. 
    :type dthresh: :py:class:`float`, optional

    :return: A filtered version of the input `reports` dictionary based 
        on the specified criteria. 
    :rtype: :py:class:`dict`
    """
    notin=[]
    for eventid in reports:
            report = reports[eventid]
            llin=[]
            llout=[]
            out=[]
            for i,lo in enumerate(report['Lon.']):
                la = report['Lat.'][i]
                if [la,lo] in llin:
                    continue
                elif [la,lo] in llout: 
                    out+=[i]
                    continue
                if isincountries(la,lo,countrycodes,dthresh=dthresh):
                    llin+=[[la,lo]]
                else:
                    llout+=[[la,lo]]
                    out+=[i]
                    
            for k in report:
                reports[eventid][k] = [e for j,e in enumerate(reports[eventid][k]) if j not in out]
    
    empty=[]
    for eventid in reports:
        if  len(reports[eventid]['Lon.'])==0:
            empty+=[eventid]
    for eventid in empty:
        reports.pop(eventid)
    return reports

def getref(reports,
           url='USGS',
           kind='us',
           dt=60*4,
           dla=1.5,
           dlo=1.5,
           f=None,
           names=[],
           sep=','):
    """
    Retrieves earthquake event data based on the input EEW report time-period and area,
    and returns the corresponding earthquake catalog.
    
    :param reports: A dictionary of EEW reports as provided by 
        :func:`pyshake.catalog.eew_report2dict` or 
        :func:`pyshake.catalog.reportfilter`.
    :type reports: :py:class:`dict`
    :param url: The source FDSN web service URL for retrieving reference earthquake 
        reports (see :external:py:mod:`obspy.clients.fdsn`).
    :type url: :py:class:`str`, optional
    :param kind: The type of web service to use for retrieving earthquake data. The 
        following are implemented: sc (SeisComP implementation of FDSN web service), or 
        us (USGS implementation of FDSN web service)
    :type kind: :py:class:`str`, optional
    :param dt: The time window in seconds that is added before and after the event 
        origin times for querying seismic events FDSN web service.
    :type dt: :py:class:`float`, optional
    :param dla: The latitude range in degrees that is added and subtracted to the event 
        latitudes for querying seismic events FDSN web service.
    :type dla: :py:class:`float`, optional
    :param dlo: The latitude range in degrees that is added and subtracted to the event 
        latitudes for querying seismic events FDSN web service.
    :type dlo: :py:class:`float`, optional
    
    :return: A reference catalog of seismic events covering the time-period and the area 
        of the input catalog.
    :rtype: :external:py:class:`pandas.DataFrame`
    """
    starttime = UTCDateTime()
    endtime = UTCDateTime('1970-01-01')
    minlon = 360
    maxlon = -360
    minlat = 90
    maxlat= -90
    for eventid in reports:
        report = reports[eventid]
        if not len(report['Lon.']):
            continue
        starttime = min([ starttime, min(report['origin time (UTC)'])])        
        endtime = max([endtime, max(report['origin time (UTC)'])])
        minlon = min([ minlon, min(report['Lon.'])])
        maxlon = max([ maxlon, max(report['Lon.'])])
        minlat = min([ minlat, min(report['Lat.'])])
        maxlat = max([ maxlat, max(report['Lat.'])])
    
    if f is None:
        f=io.BytesIO()
        opt={'starttime':starttime-dt,
            'endtime':endtime+dt,
            'minlongitude':minlon-dlo,
            'maxlongitude':maxlon+dlo,
            'minlatitude':minlat-dla,
            'maxlatitude':maxlat+dla,
            'format':'csv',
            #'limit':10,
            'filename':f}

        if kind=='us': # USGS web service
            names=['time','latitude','longitude','depth','mag',
                'magType','nst','gap','dmin','rms','net','id','updated','place','type',
                'horizontalError','depthError','magError','magNst','status',
                'locationSource','magSource']
            sep = ','
        elif kind == 'sc': # SeisComP web service
            names=['id','time','latitude','longitude','depth',
                'locationSource','Catalog','Contributor','ContributorID',
                'MagType','mag','magSource',
                'place']
            opt['format']='text'
            sep = '|'
        print(opt)
        Client(url).get_events(**opt)
    
        f.seek(0)
    else:
        print(f'Using existing {f}')

    cat = read_csv(f,header=0,names=names,sep=sep, index_col=False)
    
    if f is not None:
        cat.time=[UTCDateTime(e)+30 for e in cat.time.values]
        cat = cat[cat['time']>=starttime-dt]
        cat = cat[cat['time']<=endtime+dt]
        cat = cat[cat['longitude']>=minlon-dlo]
        cat = cat[cat['longitude']<=maxlon+dlo]
        cat = cat[cat['latitude']>=minlat-dla]
        cat = cat[cat['latitude']<=maxlat+dlo]
    else:
        if kind=='sc': # SeisComP web service
            cat = cat[cat['locationSource'].isnull() == False]
            cat = cat[cat['mag'].isnull() == False]
            cat = cat[cat['locationSource'].str.contains('scautoloc') == False]
            cat = cat[cat['locationSource'].str.contains('screloc') == False]
            cat = cat[cat['locationSource'].str.contains('scanloc') == False]
            cat = cat[cat['locationSource'].str.contains('scvsloc') == False]
            cat = cat[cat['locationSource'].str.contains('scvsnloc') == False]
        cat.time=[UTCDateTime(e) for e in cat.time.values]

    return cat
    

def matchreports(reports, cat,
                 maxhypocentraldist=150,
                 maxtimediff=100,
                 v=(5.5+3.3)/2):
    """
    Compares earthquake reports with a reference catalog using event sources parameters
    to find matching events within the specified maximum hypocentral distance and 
    origin time difference.
    
    :param reports: A dictionary of EEW reports as provided by 
        :func:`pyshake.catalog.eew_report2dict` or 
        :func:`pyshake.catalog.reportfilter`.
    :type reports: :py:class:`dict`
    :param cat: The reference catalog of seismic event as provided by 
        :func:`pyshake.catalog.getref`
    :type: :external:py:class:`pandas.DataFrame`
    :param maxhypocentraldist: The maximum hypocentral distance in kilometers allowed 
        between two seismic events to be considered as possibly matching. 
    :type maxhypocentraldist: :py:class:`float`, optional
    :param maxtimediff: The maximum time difference in seconds allowed between the origin 
        times of two seismic events to be considered as possibly matching. 
    :type maxtimediff: :py:class:`float`, optional
    :param v: The average seismic velocity used in the calculation of reduced event distance
        based on their origin time difference.
    :type v: :py:class:`float`, optional

    :return: The dictionary of EEW reports updated with the reference source parameters from 
        the matching reference seismic events.
    :rtype: :py:class:`dict`
    """
            
    if True:
        print('Similar origin times')
        ots={str(ot):[ot] for eventid in reports for ot in reports[eventid]['origin time (UTC)']}
        n=0
        for ot in ots:
            ots[ot] = cat.loc[ abs(ots[ot]-cat.time)<maxtimediff ] # () & ((ots[ot]-cat.time)>-10)
            if  len(ots[ot].longitude): 
                n+=1
        print(n,'reports with matching origin time')
    else:
        print('Similar creation times')
        cts={str(ct):[ct] for eventid in reports for ct in reports[eventid]['Creation t.']}
        n=0
        for ct in cts:
            cts[ct] = cat.loc[ ( (cts[ct]-cat.time)<maxtimediff ) & ((cat.time-cts[ct])>5) ]
            if  len(cts[ct].longitude): 
                n+=1
        print(n,'reports with matching creation time')

    print('Similar origin locations')
    n1=0
    n2=0
    for eventid in reports:
        reports[eventid]['reference']=[None for e in reports[eventid]['Lon.']]
        reports[eventid]['reference solution']=[None for e in reports[eventid]['Lon.']]

        for i,ot in enumerate(reports[eventid]['origin time (UTC)']):
            eew=[reports[eventid]['Lat.'][i],reports[eventid]['Lon.'][i],reports[eventid]['Depth'][i]]
            tmp=ots[str(ot)] 
        
        #for i,ct in enumerate(reports[eventid]['Creation t.']):
        #    ot=reports[eventid]['origin time (UTC)'][i]
        #    eew=[reports[eventid]['Lat.'][i],reports[eventid]['Lon.'][i],reports[eventid]['Depth'][i]]
        #    tmp=cts[str(ct)]

            if not len(tmp.longitude):
                continue
            hyperdistances=[]
            for j in range(len(tmp.longitude)):
                #if reports[eventid]['Creation t.'][i] < tmp.time.values[j]-5:
                #    print('EEW before origin:', reports[eventid]['Creation t.'][i], '<', tmp.time.values[j], '- 5s')
                #    continue
                ref=[tmp.latitude.values[j],tmp.longitude.values[j],tmp.depth.values[j]]
                hypocentral_distance = eventdistance(*ref,*eew)
                if hypocentral_distance>maxhypocentraldist:
                    n1+=1
                    continue
                reduced_distance = abs(tmp.time.values[j]-ot)*v
                hyper_distance = (hypocentral_distance**2+reduced_distance**2)**.5
                if False:#hyper_distance>maxhyperdist:
                    n2+=1
                    continue
                hyperdistances += [hyper_distance]
            if not len(hyperdistances):
                if reports[eventid]['Mag.'][i]>5:
                    print('Warning: missing ref for')
                    print(eventid)
                continue
            reports[eventid]['reference'][i] = tmp.id.values[argmin(hyperdistances)]
            reports[eventid]['reference solution'][i] = tmp.iloc[argmin(hyperdistances)]
    print(n1,'origin loc mismatchs')
    #print(n2,'origin loc+time mismatchs')
    return reports

def alert_accuracy(reports, references,
                 MMIcountry='Costa Rica',
                 Mtypes='Mfd,MVS',    
                 Mmin=5.5,    
                 Lmin=0.8,    
                 Dmin=0.1,    
                 MMImin=4,
                 ipm=gmm,
                 stationmin=2):
    """
    Evaluates the accuracy of earthquake reports by comparing their alert
    parameters (from the event first solution over alerting criteria) with the
    corresponding reference values based on the input reference catalog of 
    seismic event source parameters.

    The alert parameters (maximal shaking intensity, MMI) are inferred from the  
    source parameters and the alerting parameters (Mmin, Lmin, MMImin, MMIcountry, 
    ipm, and Mtypes) provided as options.
    
    :param reports: A dictionary of EEW reports with reference solutions as 
        provided by :func:`pyshake.catalog.matchreports`.
    :type reports: :py:class:`dict`
    :param references: The reference catalog of seismic events as provided by 
        :func:`pyshake.catalog.getref` to find events missing from reports.
    :type references: :external:py:class:`pandas.DataFrame`
    :param MMIcountry: The country for which the maximal Modified Mercalli Intensity 
        (MMI) is calculated. 
    :type MMIcountry: :py:class:`str`
    :param Mtypes: The comma-separated list of types of magnitudes that are 
        considered for alerting in the event reports. 
    :type Mtypes: :py:class:`str`
    :param Mmin: The minimum magnitude threshold for an earthquake event to be 
        considered in the analysis. 
    :type Mmin: :py:class:`float`, optional
    :param Lmin: The minimum likelihood value for a seismic event reports to be 
        considered in the analysis. 
    :type Lmin: :py:class:`float`, optional
    :param stationmin: The minimum number of station magnitude contributions for 
        a seismic event reports to be considered in the analysis. 
    :type stationmin: :py:class:`int`, optional
    :param Dmin: The minimum allowable difference in longitude and latitude values 
        when comparing locations. 
    :type Dmin: :py:class:`float`, optional
    :param MMImin: The minimum value of the Modified Mercalli Intensity (MMI) that 
        an earthquake must reach to be considered in the evaluation. 
    :type MMImin: :py:class:`float`, optional
    :param ipm: A prediction model to calculate the shaking intensity of a seismic event 
        based on the country distance and magnitude of the event. See :mod:`pyshake.gm`.
    :type ipm: :class:`pyshake.gm.gmm`, optional

    :return: Four dictionaries with magnitude type keys include EEW and reference 
        source parameters of Tp (true positives), Fp (false positives), Fn (false 
        negatives), and bounds.
    :rtype: :py:class:`tuple`
    """

    bounds = [180,-180,90,-90]
    Fp={mtype:[] for mtype in Mtypes.split(',')} 
    Fn={'None':[] }
    Tp={mtype:[] for mtype in Mtypes.split(',')}

    FPids=[]
    for eventid in reports:
        FP=True
        for ref in reports[eventid]['reference']:
            if ref is not None:
                FP=False
                break
        if FP and eventid not in FPids :
            FPids+=[eventid]
    FNids=[]
    for ref_id in references.id.values:
        FN=True
        for eventid in reports:
            for ref in reports[eventid]['reference']:
                if ref == ref_id:
                    FN=False
                    break
            if not FN:
                break
        if FN and ref_id not in FNids :
            FNids+=[ref_id]
    print('%d/%d F+ and %d/%d F-'%(len(FPids),len(reports),len(FNids),len(reports)))

    for evtid in reports:

        for i,mtype in enumerate(reports[evtid]['Type']):
            if Mtypes is not None and mtype not in Mtypes:
                continue
            
            nmag = reports[evtid]['Ma.'][i]
            if nmag is not None and nmag != 0 and nmag<stationmin:
                print('!!!!!',nmag,'station only in',evtid)
                continue
            norg = reports[evtid]['Or.'][i]
            if norg is not None and norg != 0 and norg<stationmin:
                print('!!!!!',norg,'station only in',evtid)
                continue

            m = reports[evtid]['Mag.'][i]
            if m is None :
                continue
            if Mmin is not None and m < Mmin:
                continue

            if reports[evtid]['Lik.'][i] is None :
                continue
            if Lmin is not None and reports[evtid]['Lik.'][i] < Lmin:
                continue
            
            lo = reports[evtid]['Lon.'][i]
            la = reports[evtid]['Lat.'][i]
            de = reports[evtid]['Depth'][i]
            di = eventcountrydistance(MMIcountry, 
                                      (lo,la,de))
            mmi = ipm.get_intensity(di,m)
            if MMImin is not None and mmi < MMImin:
                continue

            dt = reports[evtid]['Tdiff'][i]
            ot = reports[evtid]['origin time (UTC)'][i]
            eew = [lo,la,de,m,ot,dt,evtid]

            refsol = reports[evtid]['reference solution'][i]
            if evtid in FPids :
                Fp[mtype] += [eew]

            elif refsol is None or len(refsol)==0:
                print('Extra False positive: There is no reference solution for the following EEW')
                print(','.join(['%s: %s'%(k,reports[evtid][k][i]) for k in reports[evtid]]))
                Fp[mtype] += [eew]

            else:
                lo = refsol.longitude
                la = refsol.latitude
                de = refsol.depth
                m = refsol.mag
                ot = refsol.time
                dt = reports[evtid]['Creation t.'][i] - refsol.time 
                Tp[mtype] += [[lo,la,de,m,ot,dt]+eew]         

            bounds = [min([bounds[0],lo]),
                      max([bounds[1],lo]),
                      min([bounds[2],la]),
                      max([bounds[3],la])]
            break

    milo = bounds[0]
    malo = bounds[1]
    mila = bounds[2]
    mala = bounds[3]
    for ref in range(len(references.longitude)):

        if references.id.values[ref] not in FNids :
            continue

        m = references.mag.values[ref]     
        if Mmin is not None and m < Mmin:
            continue
                
        lo = references.longitude.values[ref]
        la = references.latitude.values[ref]
        if (lo<(milo-Dmin) 
            or lo>(malo+Dmin) 
            or la<(mila-Dmin) 
            or la>(mala+Dmin)):
            continue
        
        de = references.depth.values[ref]
        di = eventcountrydistance(MMIcountry, 
                                    (lo,la,de))
        mmi = ipm.get_intensity(di,m)
        if MMImin is not None and mmi < MMImin:
            continue

        ot = references.time.values[ref] 
        Fn['None'] += [[lo,la,de,m,ot,None]]
        bounds = [min([bounds[0],lo]),
                    max([bounds[1],lo]),
                    min([bounds[2],la]),
                    max([bounds[3],la])]

    return Tp, Fp, Fn, bounds
                                      
def map(reports,references,
        TFpn=None,
        rightside=True,
        top=False,
        title=None,
        fig=None,ax=None,
        figsize=(8,8),
        legendsize=True,
        legendloc=None,
        legendfaults=True,
        subtitle=False,
        colorbar=True,
        bg=True,
        **options):    
    """
    Run :func:`pyshake.catalog.alert_accuracy` (supporting all its related 
    parameters) and generate the map of results, including legends and title
    with specified options.
    
    :param reports: A dictionary of EEW reports with reference solutions as 
        provided by :func:`pyshake.catalog.matchreports`.
    :type reports: :py:class:`dict`
    :param references: The reference catalog of seismic events as provided by 
        :func:`pyshake.catalog.getref` to find events missing from reports.
    :type references: :external:py:class:`pandas.DataFrame`
    :param TFpn: Pre-computer event classification as provided by 
        :func:`pyshake.catalog.alert_accuracy`, making table output faster.
    :type TFpn: :py:class:`tuple`, optional
    :param rightside: Determines whether the legend should be displayed on the 
        right side of the plot or not (left side).
    :type rightside: :py:class:`bool`, optional
    :param title: Provide a descriptive title for the map. 
    :type title: :py:class:`str`, optional
    :param fig: Specify the figure object where the map will be plotted.
    :type fig: :external:py:class:`matplotlib.figure.Figure`, optional
    :param ax: Specify the axe object where the map will be plotted.
    :type ax: :external:py:class:`matplotlib.axes.Axes`, optional
    :param figsize: Specify the size of the figure (width, height) that will be 
        created for the plot. See :external:py:class:`matplotlib.figure.Figure`. 
    :type figsize: :py:class:`tuple`, optional
    :param legendsize: Whether the legend should include the event size scale.  
    :type legendsize: :py:class:`bool`, optional
    :param legendfaults: Whether the legend should include the fault types.  
    :type legendfaults: :py:class:`bool`, optional

    :return: The matplotlib axis where the plots have been created
    :rtype: :external:py:class:`matplotlib.axes.Axes`
    """

    if 'Mmin' not in options:
        options['Mmin']=None
    if 'Lmin' not in options:
        options['Lmin']=None
    if 'Mtypes' not in options:
        options['Mtypes']=None
    if 'MMIcountry' not in options:
        options['MMIcountry']=None
    if 'MMImin' not in options:
        options['MMImin']=None
    if 'stationmin' not in options:
        options['stationmin']=None

    if TFpn is None:
        Tp, Fp, Fn, bounds = alert_accuracy(reports,references,
                                        Mmin=options['Mmin'],
                                        Lmin=options['Lmin'],
                                        Mtypes=options['Mtypes'],
                                        MMIcountry=options['MMIcountry'],
                                        MMImin=options['MMImin'],
                                        stationmin=options['stationmin'])
    else:
        Tp, Fp, Fn, bounds = TFpn

    n = 0
    for tmp in [Tp, Fn]:
        for mtype in tmp:
            n+=len(tmp[mtype])
    
    ax=bmap(bounds=bounds,
            fig=fig,
            ax=ax,
            figsize=figsize,
            legendfaults=legendfaults,
            bg=bg,
            top=top,
            rightside=rightside)

    times = [xyz[4] for mtype in Tp for xyz in Tp[mtype]]+\
            [xyz[4] for mtype in Fp for xyz in Fp[mtype]]+\
            [xyz[4] for mtype in Fn for xyz in Fn[mtype]]
    
    times = [min(times),max(times)]

    mtypes = list(set([mtype for mtype in Tp]+[mtype for mtype in Fp])) + \
             ['None']
    
    scatter,scattercmaps = plot(Tp,
                   mtypes,
                   ax,
                   marker='o',
                   times=times,
                   transform=ccrs.Geodetic(),
                   zorder=97)
    
    plot({mtype: [ xyz[:6] for xyz in Fp[mtype]] for mtype in Fp},
         mtypes,
         ax,
         label='Fp$^{%s}$',
         marker='X',
         times=times,
         transform=ccrs.Geodetic(),
         nodt=True,
         zorder=98)
    
    plot(Fn,
         mtypes,
         ax,
         label='Fn$^{%s}$',
         marker='s',
         times=times,
         transform=ccrs.Geodetic(),
         zorder=99)
    
    ax.legend()
    h, l = ax.get_legend_handles_labels()

    if legendsize:
        kw = dict(prop="sizes", 
            num=4, 
            color='k',#scatter[0].cmap(scatter[0].norm(150)), 
            fmt="M{x:.2g}",
            func=size2mag)
        h2, l2 = [l[::-1] for l in scatter[0].legend_elements(**kw)]

        h = h2 + h[::-1]
        l = l2 + l[::-1]
        h = h[::-1]
        l = l[::-1]

    locbar = 'left'
    loclegend = 'right'
    bbox_to_anchor = 0
    
        
    if colorbar:    

        bbox = ax.get_position() 

        if not rightside:
            locbar = 'right'
            loclegend = 'left'
            bbox_to_anchor = 1
            cax = ax.figure.add_axes([bbox.x0+bbox.width, 
                                    bbox.y0+bbox.height-bbox.height/2.5, 
                                    bbox.width/1.5, 
                                    bbox.height/2.5])
            
        else:
            cax = ax.figure.add_axes([bbox.x0-bbox.width/2.5,#-bbox.width/1.5, 
                                    bbox.y0+bbox.height/15,#.5bbox.height-bbox.height/2.5, 
                                    bbox.width/1.5, 
                                    bbox.height/3])
            
            #cax = inset_axes(ax,
            #                 width="5%",  # width: 50% of parent_bbox width
            #                 height="5%",  # height: 5%
            #                 loc="lower left",
            #                 )
        cax.axis('off')

        for im in scattercmaps:
            cbar = ax.figure.colorbar(im[0], 
                            ax=cax,
                            location=locbar,
                            fraction=0.8,
                            pad=0,
                            panchor=(0,0),
                            #use_gridspec=True,
                            #shrink=0.5
                            )
            cbar.set_label(im[1], 
                            rotation=90)
            if locbar == 'left':
                locbar = 'right'
            else:
                locbar = 'left'


    l=ax.legend(h, l,
                #ncol=1,
                #fontsize='x-small',
                #bbox_to_anchor=(bbox_to_anchor,0), 
                loc=legendloc,
                )        
    
    
    t = '%d ev.'%n

    if title is not None:
        t = '%s %s'%(title,t)
    
    if subtitle:
        t += '\n$^{%s}_{%s}$'%(times[0].isoformat()[:16], times[1].isoformat()[:16])
        
        if options['MMImin'] is not None:
            t += '\n$MMI$'
            if options['MMIcountry'] is not None:
                t += '$_{%s}$'%options['MMIcountry'].capitalize() 
            t += '$>%.1g$'%options['MMImin'] 

        if options['Mmin'] is not None and options['Mtypes'] is not None:
            t += '\n$M_{%s}$'%options['Mtypes'].replace('M','')
            t += '$>%.1f$'%options['Mmin']

        if options['Lmin'] is not None:
            t += '\n$Lik.>%.1g$'%options['Lmin']
        
    
    l.set_title(title=t, prop={'weight':'bold'})

    
    return ax



def table(reports=None,references=None,TFpn=None,**options):
    """
    Run :func:`pyshake.catalog.alert_accuracy` (supporting all its related 
    parameters) and generate a table of false events.
    
    :param reports: A dictionary of EEW reports with reference solutions as 
        provided by :func:`pyshake.catalog.matchreports`.
    :type reports: :py:class:`dict`
    :param references: The reference catalog of seismic events as provided by 
        :func:`pyshake.catalog.getref` to find events missing from reports.
    :type references: :external:py:class:`pandas.DataFrame`
    :param TFpn: Pre-computer event classification as provided by 
        :func:`pyshake.catalog.alert_accuracy`, making table output faster.
    :type TFpn: :py:class:`tuple`, optional

    :return: The matplotlib axis where the plots have been created
    :rtype: :external:py:class:`matplotlib.axes.Axes`
    """
    
    if 'Mmin' not in options:
        options['Mmin']=None
    if 'Lmin' not in options:
        options['Lmin']=None
    if 'Mtypes' not in options:
        options['Mtypes']=None
    if 'MMIcountry' not in options:
        options['MMIcountry']=None
    if 'MMImin' not in options:
        options['MMImin']=None
    if 'stationmin' not in options:
        options['stationmin']=None

    if TFpn is None:
        Tp, Fp, Fn, bounds = alert_accuracy(reports,references,
                                        Mmin=options['Mmin'],
                                        Lmin=options['Lmin'],
                                        Mtypes=options['Mtypes'],
                                        MMIcountry=options['MMIcountry'],
                                        MMImin=options['MMImin'],
                                        stationmin=options['stationmin'])
    else:
         Tp, Fp, Fn, bounds = TFpn
    #print(Fp)
    #print(Fn)
    
    columns=('Longitude','Latitude','Depth','Magnitude','Origin time')#,'EEW delay')
    dtfs = []
    for mtype in Fp:
        dtfs += [DataFrame([xyz[:5] for xyz in Fp[mtype]], #Fp[mtype],
                            columns=columns,
                            index=[('F+',mtype) for e in Fp[mtype] ])]
    
    for mtype in Fn:
        dtfs += [DataFrame([xyz[:5] for xyz in Fn[mtype]], #Fn[mtype],
                            columns=columns,
                            index=[('F-','') for e in Fn[mtype] ])]
        
    return concat(dtfs)


def errors(reports=None,
           references=None,
           TFpn=None,
           ax=None,
           title=None,
           subtitle=False,
           cmaps={'Mfd':'autumn','MVS':'winter'},
           colorbar=True,
           rightside=False,
           cax=None,
           legendsize=True,
           lax=None,
           **options):
    """
    Run :func:`pyshake.catalog.alert_accuracy` (supporting all its related 
    parameters) and generate a location and magnitude error histogram.
    
    :param reports: A dictionary of EEW reports with reference solutions as 
        provided by :func:`pyshake.catalog.matchreports`.
    :type reports: :py:class:`dict`
    :param references: The reference catalog of seismic events as provided by 
        :func:`pyshake.catalog.getref` to find events missing from reports.
    :type references: :external:py:class:`pandas.DataFrame`
    :param TFpn: Pre-computer event classification as provided by 
        :func:`pyshake.catalog.alert_accuracy`, making table output faster.
    :type TFpn: :py:class:`tuple`, optional

    :return: The table of false events that can be viewed with :py:func:`print` 
        or :py:func:`display`.
    :rtype: :external:py:class:`matplotlib.DataFrame`
    """
    
    if ax is None:
        ax = figure().gca()
        
    if 'Mmin' not in options:
        options['Mmin']=None
    if 'Lmin' not in options:
        options['Lmin']=None
    if 'Mtypes' not in options:
        options['Mtypes']=None
    if 'MMIcountry' not in options:
        options['MMIcountry']=None
    if 'MMImin' not in options:
        options['MMImin']=None
    if 'stationmin' not in options:
        options['stationmin']=None

    if TFpn is None:
        Tp, Fp, Fn, bounds = alert_accuracy(reports,references,
                                        Mmin=options['Mmin'],
                                        Lmin=options['Lmin'],
                                        Mtypes=options['Mtypes'],
                                        MMIcountry=options['MMIcountry'],
                                        MMImin=options['MMImin'],
                                        stationmin=options['stationmin'])
    else:
         Tp, Fp, Fn, bounds = TFpn

    def _forward(x):
        return sqrt(x)
    def _inverse(x):
        return x**2

    norm = colors.FuncNorm((_forward, _inverse), 
                            vmin=0, 
                            vmax=200)
    
    n = 0 
    scattercmaps={}
    beginend = sort([ xyzmte[4] for mtype in Tp for xyzmte in Tp[mtype]  ])
    for i,mtype in enumerate(Tp):
        
        n += len(Tp[mtype])
        
        dist = [ eventcountrydistance(options['MMIcountry'], (xyzmte[0],xyzmte[1],xyzmte[2]) ) for xyzmte in Tp[mtype]  ]
        d_epi = [ gps2dist_azimuth(xyzmte[1],xyzmte[0],xyzmte[7],xyzmte[6])[0]/1000.0 for xyzmte in Tp[mtype]  ] 
        d_dep = [ xyzmte[8] - xyzmte[2] for xyzmte in Tp[mtype]  ] 
        d_hyp = [(d_epi[i]**2 + d_dep[i]**2)**.5 for i,x in enumerate(d_dep)]
        d_mag = [ xyzmte[9]-xyzmte[3] for xyzmte in Tp[mtype] ] 
        mag = [ xyzmte[3] for xyzmte in Tp[mtype] ]     

        times = [ xyzmte[4] for xyzmte in Tp[mtype] ]   
        alphas = asarray(sortarg(times)) 
        alphas = alphas + max(alphas)/3 
        alphas = alphas / max(alphas)
        
        c = cm.ScalarMappable(norm=norm, cmap=cmaps[mtype]).to_rgba(dist)#[xyzmte[2] for xyzmte in Tp[mtype]])#
        c[:,-1] = alphas
        
        scattercmaps[mtype] = ax.scatter([nan,nan],[nan,nan],
                                         mag2size(asarray([3,7])),
                                         [nan,nan],
                                         cmap=cmaps[mtype],   
                                         norm=norm,
                                         alpha=0.5)
        
        ax.scatter(d_hyp,
                    d_mag,
                    mag2size(asarray(mag)),
                    c,#dist,
                    marker='oo'[i],           
                    linewidths=1.0,
                    zorder=9999)

    locbar = 'left' 
    bbox = ax.get_position() 
    if cax is None:
        if rightside:
            locbar = 'right'
            cax = ax.figure.add_axes([bbox.x0+bbox.width, 
                                    bbox.y0+bbox.height-bbox.height/2.5, 
                                    bbox.width/1.5, 
                                    bbox.height/2.5])
        else:
            cax = ax.figure.add_axes([bbox.x0-bbox.width/1,#.5, 
                                    bbox.y0+bbox.height-bbox.height/2.5, 
                                    bbox.width/1.5, 
                                    bbox.height/2.5])
    cax.axis('off')
        
    if len(scattercmaps)==1:
        locbar = 'right'

    if colorbar:    
        for mtype in scattercmaps:
            cbar = ax.figure.colorbar(scattercmaps[mtype], 
                            ax=cax,
                            location=locbar,
                            fraction=.8,
                            pad=0,
                            panchor=(0,0),
                            #use_gridspec=True,
                            #shrink=0.5
                            )
            cbar.set_label('M$_{%s}$'%mtype[1:], 
                            rotation=90)
            if locbar == 'left':
                locbar = 'right'
            else:
                locbar = 'left'
            cbar.solids.set(alpha=1)

        cax.set_title('Country\ndistance (km)',#Depth (km)',#
                      loc=['left','center'][len(Tp)>1])

    l=ax.legend()        

    t = '%d Tp'%( n )

    if title is not None:
        t = '%s%s-%s\n%s'%(title, 
                           str(beginend[0])[:7].replace('-','/'), 
                           str(beginend[-1])[:7].replace('-','/'), 
                           t)
    
    if subtitle and 'MMImin' in options :
        t += '\n$^{%s}_{%s}$'%(sorted(times)[0].isoformat()[:16], sorted(times)[1].isoformat()[:16])
        
        if options['MMImin'] is not None:
            t += '\n$MMI$'
            if options['MMIcountry'] is not None:
                t += '$_{%s}$'%options['MMIcountry'].capitalize() 
            t += '$>%.1g$'%options['MMImin'] 

        if options['Mmin'] is not None and options['Mtypes'] is not None:
            t += '\n$M_{%s}$'%options['Mtypes'].replace('M','')
            t += '$>%.1f$'%options['Mmin']

        if options['Lmin'] is not None:
            t += '\n$Lik.>%.1g$'%options['Lmin']        
    
    l.set_title(title=t, prop={'weight':'bold'})

    if legendsize:
        ax.add_artist(l)
        kw = dict(prop="sizes", 
                    num=[4,5,6], 
                    alpha=0.8,
                    color='k',
                    fmt="M{x:.2g}",
                    func=size2mag)
        if lax is None:
            lax = ax
        l = lax.legend(*scattercmaps[list(scattercmaps.keys())[0]].legend_elements(**kw),
                            #loc="lower right", 
                            title="Magnitude",
                            loc='lower left')
        lax.add_artist(l)
        lax.scatter([nan],[nan],
                    mag2size(asarray([5])),
                    'k',
                    label='Earliest',
                    alpha=0.3)
        lax.scatter([nan],[nan],
                    mag2size(asarray([5])),
                    'k',
                    label='Latest',
                    alpha=1)
        l = lax.legend(title="Time",
                            loc='upper left')
        lax.add_artist(l)


    ax.tick_params(right=True, top=True,
                    left=True, bottom=True,
                    which='both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid( which='major', color='gray', linestyle='dashdot', zorder=-9999) #b=True,
    ax.grid( which='minor', color='beige',  ls='-', zorder=-9999) #b=True,

    return ax

def event_eew_mag(events,
                  MMIcountry='Costa Rica',
                  Mtypes='Mfd,MVS',    
                  Mmin=5.5,    
                  Lmin=0.8,      
                  MMImin=4,
                  ipm=gmm,
                  debug=False):
    """
    Evaluates alert parameters (from the event first solution over alerting criteria) 
    on the input catalog of seismic event source parameters.

    The alert parameters (maximal shaking intensity, MMI) are inferred from the  
    source parameters and the alerting parameters (Mmin, Lmin, MMImin, MMIcountry, 
    ipm, and Mtypes) provided as options.
    
    :param events: The reference catalog of seismic events.
    :type events: :external:py:class:`obspy.core.event.Catalog`
    :param MMIcountry: The country for which the maximal Modified Mercalli Intensity 
        (MMI) is calculated. 
    :type MMIcountry: :py:class:`str`
    :param Mtypes: The comma-separated list of types of magnitudes that are 
        considered for alerting in the event reports. 
    :type Mtypes: :py:class:`str`
    :param Mmin: The minimum magnitude threshold for an earthquake event to be 
        considered in the analysis. 
    :type Mmin: :py:class:`float`, optional
    :param Lmin: The minimum likelihood value for a seismic event reports to be 
        considered in the analysis. 
    :type Lmin: :py:class:`float`, optional
    :param Dmin: The minimum allowable difference in longitude and latitude values 
        when comparing locations. 
    :type Dmin: :py:class:`float`, optional
    :param MMImin: The minimum value of the Modified Mercalli Intensity (MMI) that 
        an earthquake must reach to be considered in the evaluation. 
    :type MMImin: :py:class:`float`, optional
    :param ipm: A prediction model to calculate the shaking intensity of a seismic event 
        based on the country distance and magnitude of the event. See :mod:`pyshake.gm`.
    :type ipm: :class:`pyshake.gm.gmm`, optional

    :return: First event magnitudes exceeding the alert thresholds.
    :rtype: :py:class:`list`
    """
    
    eew_mags = [ None for e in events]
    
    for evtindex,e in enumerate(events):
        
        if e is None:
            continue

        ct = [m.creation_info.creation_time for m in  e.magnitudes]
        for m in argsort(ct):

            if eew_mags[evtindex] is not None:
                print('This should never append')
                continue

            magnitude = e.magnitudes[m]
            origin = ([None]+[o for o in e.origins if magnitude.origin_id.id == o.resource_id.id])[-1]

            if debug:
                print(origin.longitude,origin.latitude,origin.depth/1000.,magnitude.mag,origin.time)

            if origin is None:
                continue
            
            if Mtypes is not None and magnitude.magnitude_type not in Mtypes:
                if debug:
                    print(magnitude.magnitude_type,'not in',Mtypes)
                continue

            if magnitude.mag is None :
                if debug:
                    print(magnitude.mag,'is', None )
                continue

            if Mmin is not None and magnitude.mag < Mmin:
                if debug:
                    print(magnitude.mag,'<', Mmin)
                continue

            likelihood = [None]
            likelihood += [float(c.text) for c in magnitude.comments if c.resource_id.id.split('/')[-1] == 'likelihood']
            likelihood = likelihood[-1]
            if likelihood is None :
                if debug:
                    print(likelihood,'is', None)
                continue

            if Lmin is not None and likelihood < Lmin:
                if debug:
                    print(likelihood,'<', Lmin)
                continue
            
            di = eventcountrydistance(MMIcountry, 
                                      (origin.longitude,
                                       origin.latitude,
                                       origin.depth/1000.))
            mmi = ipm.get_intensity(di,
                                    magnitude.mag)
            if MMImin is not None and mmi < MMImin:
                if debug:
                    print(mmi,' <', MMImin)
                continue
            
            eew_mags[evtindex] = [magnitude, origin]
            break

    return eew_mags

def chronology(reports=None,
              references=None,
              TFpn=None,
              ax=None,
              title=None,
              subtitle=False,
              cmaps={'Mfd':'autumn','MVS':'winter','None':None},
              colorbar=True,
              rightside=False,
              cax=None,
              legendsize=True,
              lax=None,
              xml=True,
              arrivalrank=4,
              arrivalF=False,
              legend_marker_size = 70,
              legend_marker_alpha = 2/3,
              hlegend_line_length = 1,
              annotation = None,
              **options):
    """
    Run :func:`pyshake.catalog.alert_accuracy` (supporting all its related 
    parameters) and generate a chronology plot.
    
    :param reports: A dictionary of EEW reports with reference solutions as 
        provided by :func:`pyshake.catalog.matchreports`.
    :type reports: :py:class:`dict`
    :param references: The reference catalog of seismic events as provided by 
        :func:`pyshake.catalog.getref` to find events missing from reports.
    :type references: :external:py:class:`pandas.DataFrame`
    :param TFpn: Pre-computer event classification as provided by 
        :func:`pyshake.catalog.alert_accuracy`, making table output faster.
    :type TFpn: :py:class:`tuple`, optional

    :return: The table of false events that can be viewed with :py:func:`print` 
        or :py:func:`display`.
    :rtype: :external:py:class:`matplotlib.DataFrame`
    """
    
    if ax is None:
        ax = figure().gca()
        
    if 'Mmin' not in options:
        options['Mmin']=None
    if 'Lmin' not in options:
        options['Lmin']=None
    if 'Mtypes' not in options:
        options['Mtypes']=None
    if 'MMIcountry' not in options:
        options['MMIcountry']=None
    if 'MMImin' not in options:
        options['MMImin']=None
    if 'stationmin' not in options:
        options['stationmin']=None

    if TFpn is None:
        Tp, Fp, Fn, bounds = alert_accuracy(reports,references,
                                        Mmin=options['Mmin'],
                                        Lmin=options['Lmin'],
                                        Mtypes=options['Mtypes'],
                                        MMIcountry=options['MMIcountry'],
                                        MMImin=options['MMImin'],
                                        stationmin=options['stationmin'])
    else:
         Tp, Fp, Fn, bounds = TFpn

    def _forward(x):
        return sqrt(x)
    def _inverse(x):
        return x**2

    norm = colors.FuncNorm((_forward, _inverse), 
                            vmin=0, 
                            vmax=200)
    
    def linnorm(x):
        y = 1 - asarray(x)/150 #asarray(sortarg(times)) + 1
        try:
            y[y>1]=1
            y[y<0]=0
        except:
            pass
        y = y * 3/4 + 1/4 #/ max(alphas)
        return y

    scattercmaps={}
    beginend = sort([ xyzmte[4] for mtype in Tp for xyzmte in Tp[mtype]  ])
    n=0
    for j,cl in  enumerate([Tp,Fp,Fn]):
        for i,mtype in enumerate(cl):
   
            dist = [ eventcountrydistance(options['MMIcountry'], (xyzmte[0],xyzmte[1],xyzmte[2]) ) for xyzmte in cl[mtype]  ]
            
            if j==0:
                d_epi = [ gps2dist_azimuth(xyzmte[1],xyzmte[0],xyzmte[7],xyzmte[6])[0]/1000.0 for xyzmte in cl[mtype]  ] 
                d_dep = [ xyzmte[8] - xyzmte[2] for xyzmte in cl[mtype]  ] 
                d_hyp = [(d_epi[i]**2 + d_dep[i]**2)**.5 for i,x in enumerate(d_dep)]
            else:
                d_hyp = [ 0 for xyzmte in cl[mtype]  ] 
            
            
            d_time = [ xyzmte[11-(j*6)] for xyzmte in cl[mtype]  ] 
            if j==2 or j==1:
                d_time = [ 0 for xyzmte in cl[mtype]  ] 

            mag = [ xyzmte[3] for xyzmte in cl[mtype] ]      
            times = [ xyzmte[4].datetime for xyzmte in cl[mtype] ]          
            arrivals = [ nan for xyzmte in cl[mtype] ]         
            
            if xml and ( arrivalF or j == 0 ):
                
                for e,xyzmte in enumerate(cl[mtype]):  
                    
                    if xyzmte[14-(j*6)] is None:
                        print(14-(j*6),xyzmte)
                        continue


                    magnitude = xyzmte[14-(j*6)][0]

                    test = magnitude.creation_info.creation_time - xyzmte[10-(j*6)]                    
                    d_time[e] = test
                    
                    event = xyzmte[15-(j*6)][0]
                    if len(event.origins)>1:
                        print('Not sure to get preferred origin for ')
                        print(event)
                    preforgtime = event.origins[0].time

                    arrivaltimes = {}
                    
                    for p in event.picks:
                        station = '%s.%s'%(p.waveform_id['network_code'],p.waveform_id['station_code'])
                        if station not in arrivaltimes:
                            arrivaltimes[station] = p.time
                        if p.time < arrivaltimes[station] :
                            arrivaltimes[station] = p.time

                    test = [arrivaltimes[s] for s in arrivaltimes]

                    if len(test)>=arrivalrank:
                        arrivals[e] = sort(test)[arrivalrank-1] - preforgtime # xyzmte[10-(j*6)] 

                
            alphas = linnorm(dist)
            
            c = cm.ScalarMappable(norm=norm, cmap=cmaps[mtype]).to_rgba(d_hyp)#dist)
            c[:,-1] = alphas

            if j == 2:
                c[:,0] = 1
                c[:,1] = 0
                c[:,2] = 1

            if j==0 or j==2:
                n += len(cl[mtype])
            
            if j==0 :
                scattercmaps[mtype] = ax.scatter([nan,nan],[nan,nan],
                                                mag2size(asarray([3,7])),
                                                [nan,nan],
                                                cmap=cmaps[mtype],   
                                                norm=norm,        
                                                alpha=0.5)
                
                for i,dt in enumerate(d_time):
                    if dt < 2 :
                        d_time[i] = nan
                        print('Ignoring ref/EEW matching issue with:')
                        print(cl[mtype][i])
            
            sc = ax.scatter(times,
                    d_time,
                    mag2size(asarray(mag)),
                    c,
                    marker='oXs'[j],           
                    linewidths=1.0,
                    label=('%s (%d)'%(['Tp$^{%s}$','Fp$^{%s}$','Fn$^{%s}$'][j]%(mtype[1:]), len(cl[mtype]))).replace('$^{one}$',''),
                    clip_on = j==0,
                    zorder=9999)
            
            ax.scatter(times,
                    arrivals,
                    mag2size(asarray(mag)),
                    c,
                    marker='___'[j],           
                    linewidths=2.5,
                    zorder=9999)
        
        
    ax.scatter([nan,nan],
                [nan,nan],
                mag2size(asarray([3,7])),
                'k',
                marker='_',           
                linewidths=1.0,
                label='%d$^{th}$ arrival'%(arrivalrank),
                zorder=9999)
    
    ax.set_ylim(bottom=0)
    bbox = ax.get_position() 
    if cax is None:
        if rightside:
            locbar = 'right'
            cax = ax.figure.add_axes([bbox.x0+bbox.width, 
                                    bbox.y0+bbox.height-bbox.height/2.5, 
                                    bbox.width/1.5, 
                                    bbox.height/2.5])
        else:
            cax = ax.figure.add_axes([bbox.x0-bbox.width/1.5, 
                                    bbox.y0+bbox.height-bbox.height/2.5, 
                                    bbox.width/1.5, 
                                    bbox.height/2.5])
    cax.axis('off')
        
    locbar = 'right'
    if colorbar:    
        for mtype in scattercmaps:
            cbar = ax.figure.colorbar(scattercmaps[mtype], 
                            ax=cax,
                            location=locbar,
                            fraction=.8,
                            pad=0,
                            panchor=(0,0),
                            #use_gridspec=True,
                            #shrink=0.5
                            )
            cbar.set_label('M$_{%s}$'%mtype[1:], 
                            rotation=90)
            if locbar == 'left':
                locbar = 'right'
            else:
                locbar = 'left'
            cbar.solids.set(alpha=1)

        cax.set_title('Location\nerror (km, hyp.)',
                      loc=['left','center'][len(Tp)>1])

    if annotation is not None:
        for k in annotation:
            ax.axvline(UTCDateTime(annotation[k]),label=k, color='k', linewidth=3, alpha=0.5, zorder = -99)   

    l=ax.legend(handlelength=hlegend_line_length)        

    for lh in l.legendHandles:
        lh.set_alpha(legend_marker_alpha)
        if hasattr(lh,'set_sizes'):
            lh.set_sizes([legend_marker_size])


    t = '%d ev.'%( n )

    if title is not None:
        t = '%s%s-%s\n%s'%(title, 
                           str(beginend[0])[:7].replace('-','/'), 
                           str(beginend[-1])[:7].replace('-','/'), 
                           t)
    
    if subtitle and 'MMImin' in options :
        t += '\n$^{%s}_{%s}$'%(sorted(times)[0].isoformat()[:16], sorted(times)[1].isoformat()[:16])
        
        if options['MMImin'] is not None:
            t += '\n$MMI$'
            if options['MMIcountry'] is not None:
                t += '$_{%s}$'%options['MMIcountry'].capitalize() 
            t += '$>%.1g$'%options['MMImin'] 

        if options['Mmin'] is not None and options['Mtypes'] is not None:
            t += '\n$M_{%s}$'%options['Mtypes'].replace('M','')
            t += '$>%.1f$'%options['Mmin']

        if options['Lmin'] is not None:
            t += '\n$Lik.>%.1g$'%options['Lmin']        
    
    l.set_title(title=t, prop={'weight':'bold'})

    if legendsize:
        ax.add_artist(l)
        kw = dict(prop="sizes", 
                    num=[4,5,6], 
                    alpha=0.8,
                    color='k',
                    fmt="M{x:.2g}",
                    func=size2mag)
        if lax is None:
            lax = ax
        l = lax.legend(*scattercmaps[list(scattercmaps.keys())[0]].legend_elements(**kw),
                            #loc="lower right", 
                            title="Magnitude",
                            loc='lower left')
        lax.add_artist(l)
        
        for d in [150,75,0]:
            lax.scatter([nan],[nan],
                        mag2size(asarray([5])),
                        'k',
                        label='%d'%d,
                        alpha=linnorm(d))
            
        l = lax.legend(title="Country\ndistance\n(km)",
                            loc='upper left')
        lax.add_artist(l)

    if ax.get_ylim()[1]>100:
        ax.set_yscale('symlog',linthresh=100)
    ax.tick_params(right=True, top=True,
                    left=True, bottom=True,
                    which='both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid( which='major', color='gray', linestyle='dashdot', zorder=-9999) #b=True,
    ax.grid( which='minor', color='beige',  ls='-', zorder=-9999) #b=True,
     
    return ax