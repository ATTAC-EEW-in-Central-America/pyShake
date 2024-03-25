# -*- coding: utf-8 -*-
"""
Catalog module. 
"""

import io, glob

from numpy import argmin 
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from pandas import read_csv, DataFrame, concat
import cartopy.crs as ccrs

from pyshake.util import isincountries, eventdistance, eventcountrydistance 
from pyshake.plotting import plot, bmap, size2mag
from pyshake.gm import gmm

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
        should returned as a :external:py:mod:`pandas.Dataframe` or a
        :external:py:mod:`obspy.Catalog`.
    :type string: :py:class:`string`.

    .. Note::

        Any additional keyword arguments will be passed to the web-service as
        additional arguments. Passing any non-default parameters that the 
        web-service does not support will raise an error.

    :return: The list of events.
    :rtype: :external:py:mod:`pandas.Dataframe` or :external:py:mod:`obspy.Catalog`
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
                    if e==True:
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
           dlo=1.5):
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

    #to do :
    #- eval Mc
    
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
        maxlat= max([ maxlat, max(report['Lat.'])])

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

    if kind == 'us': # USGS web service
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
        opt['format'] = 'text'
        sep = '|'
    
    Client(url).get_events(**opt)
    
    f.seek(0)
    cat = read_csv(f,header=0,names=names,sep=sep, index_col=False)
    
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
                 maxtimediff=3*60,
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
            
    print('Similar origin times')
    ots={str(ot):[ot] for eventid in reports for ot in reports[eventid]['origin time (UTC)']}
    n=0
    for ot in ots:
        ots[ot] = cat.loc[abs(cat.time-ots[ot])<maxtimediff]
        if  len(ots[ot].longitude): 
            n+=1
    print(n,'reports with matching origin time')

    print('Similar origin locations')
    n1=0
    n2=0
    for eventid in reports:
        reports[eventid]['reference']=[None for e in reports[eventid]['Lon.']]
        reports[eventid]['reference solution']=[None for e in reports[eventid]['Lon.']]
        for i,ot in enumerate(reports[eventid]['origin time (UTC)']):
            eew=[reports[eventid]['Lat.'][i],reports[eventid]['Lon.'][i],reports[eventid]['Depth'][i]]
            tmp=ots[str(ot)]
            if not len(tmp.longitude):
                continue
            hyperdistances=[]
            for j in range(len(tmp.longitude)):
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
                 ipm=gmm):
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
            eew = [lo,la,de,m,ot,dt]

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
        title=None,
        fig=None,ax=None,
        figsize=(8,8),
        legendsize=True,
        legendfaults=True,
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

    if TFpn is None:
        Tp, Fp, Fn, bounds = alert_accuracy(reports,references,
                                        Mmin=options['Mmin'],
                                        Lmin=options['Lmin'],
                                        Mtypes=options['Mtypes'],
                                        MMIcountry=options['MMIcountry'],
                                        MMImin=options['MMImin'])
    else:
        Tp, Fp, Fn, bounds = TFpn
    n = 0
    for tmp in [Tp, Fp, Fn]:
        for mtype in tmp:
            n+=len(tmp[mtype])
    
    ax=bmap(bounds=bounds,fig=fig,ax=ax,figsize=figsize,legendfaults=legendfaults,rightside=rightside)

    times =  [xyz[4] for mtype in Tp for xyz in Tp[mtype]]+[xyz[4] for mtype in Fp for xyz in Fp[mtype]]+[xyz[4] for mtype in Fn for xyz in Fn[mtype]]
    times = [min(times),max(times)]

    mtypes =  list(set([mtype for mtype in Tp]+[mtype for mtype in Fp]))+['None']
    scatter = plot(Tp,mtypes,ax,marker='o',times=times,transform=ccrs.Geodetic(),zorder=97)
    plot(Fp,mtypes,ax,label='F$^{+}_{%s}$',marker='X',times=times,transform=ccrs.Geodetic(),zorder=98)
    plot(Fn,mtypes,ax,label='F$^{-}_{%s}$',marker='s',times=times,transform=ccrs.Geodetic(),zorder=99)
    
    ax.legend()
    h, l = ax.get_legend_handles_labels()

    if legendsize:
        kw = dict(prop="sizes", 
            num=4, 
            color='k',#scatter[0].cmap(scatter[0].norm(150)), 
            fmt="M{x:.2g}",
            func=size2mag)
        h2, l2 = [l[::-1] for l in scatter[0].legend_elements(**kw)]
        h = h[::-1]+h2
        l = l[::-1]+l2
    if rightside:
        l=ax.legend(h, l,
                bbox_to_anchor=(1,0), 
                loc="lower left")
    else:
        l=ax.legend(h, l,
                bbox_to_anchor=(0,0), 
                loc="lower right")

    t = '%d ev.'%n

    if title is not None:
        t = '%s %s'%(title,t)
    
    if options['Mmin'] is not None and options['Mtypes'] is not None:
        t += '\n$M_{%s}$'%options['Mtypes'].replace('M','')
        t += '$>%.1f$'%options['Mmin']
    
    if options['MMImin'] is not None:
        t += '\n$MMI$'
        if options['MMIcountry'] is not None:
            t += '$_{%s}$'%options['MMIcountry'].capitalize() 
        t += '$>%.1g$'%options['MMImin'] 

    if options['Lmin'] is not None:
        t += '\n$Lik.>%.1g$'%options['Lmin']
    
    t += '\n$^{%s}_{%s}$'%(times[0].isoformat()[:16], times[1].isoformat()[:16])
    
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

    :return: The table of false events that can be viewed with :py:func:`print` 
        or :py:func:`display`.
    :rtype: :external:py:class:`pandas.Dataframe`
    """
    
    if TFpn is None:
        Tp, Fp, Fn, bounds = alert_accuracy(reports,references,
                                        Mmin=options['Mmin'],
                                        Lmin=options['Lmin'],
                                        Mtypes=options['Mtypes'],
                                        MMIcountry=options['MMIcountry'],
                                        MMImin=options['MMImin'])
    else:
         Tp, Fp, Fn, bounds = TFpn
    #print(Fp)
    #print(Fn)
    
    columns=('Longitude','Latitude','Depth','Magnitude','Origin time','EEW delay')
    dtfs = []
    for mtype in Fp:
        dtfs += [DataFrame(Fp[mtype],
                            columns=columns,
                            index=['%s & F+'%mtype for e in Fp[mtype] ])]
    
    for mtype in Fn:
        dtfs += [DataFrame(Fn[mtype],
                            columns=columns,
                            index=['F-' for e in Fn[mtype] ])]
        
    return concat(dtfs)