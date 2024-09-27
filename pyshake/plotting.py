# -*- coding: utf-8 -*-
"""
Plotting utility module. 
"""

import glob, matplotlib.pyplot, matplotlib.patheffects, matplotlib.ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import cartopy.io.shapereader as shpreader
import geopandas,requests,glob,os
import matplotlib.transforms as mtransforms
from matplotlib import colors
from numpy import argsort
from numpy import sqrt


#def argsort(seq):
#    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
#   return sorted(range(len(seq)), key=seq.__getitem__)

def mag2size(mag):
    return ((mag+.5)**5.3)/7/14
def size2mag(size):
    return (size*7*14)**(1/5.3)-.5

def mapcities(bmap,
              cities={'San Jose':[9.94456,-84.11525],
                      'Managua':[12.12801,-86.29536],
                      'San Salvador':[13.70178, -89.21693],
                      'Guatemala City':[14.62778, -90.51520]},
              optcities={#'weight':"bold",
                         'color':"k",
                         'fontsize':'x-small',
                         'zorder':9999999999},
              optcitydot={'markeredgecolor':'w',
                          'markeredgewidth':.8,
                          'markerfacecolor':'none',
                          'zorder':9999999999},
                          **opt):
    dot_path_effects=[matplotlib.patheffects.withStroke(linewidth=2,foreground="k")]
    label_path_effects=[matplotlib.patheffects.withStroke(linewidth=2,foreground="w")]
    for city in cities:
        text = bmap.plot(*cities[city][::-1],'o',path_effects=dot_path_effects,**optcitydot,**opt)
        text = bmap.text(*[d-0.02 for d in cities[city][::-1]],city,va='top',ha='right',path_effects=label_path_effects,clip_on=True,**optcities,**opt)

def mapcountrynames(ax,**opt):
    extent=ax.get_extent(opt['transform'])
    shpfilename = shpreader.natural_earth(resolution='110m',
                                      category='cultural',
                                      name='admin_0_countries')

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()
    for country in countries:

        x = country.geometry.centroid.x        
        y = country.geometry.centroid.y

        if  x<extent[0] or x>extent[1]:
            continue
        if y<extent[2] or y>extent[3]:
            continue

        ax.text(x, y, country.attributes['NAME'], 
                color='k', ha='center', va='bottom', alpha=.7,clip_on=True,
                path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w", alpha=.7)], 
                **opt)

def sanitize_lonlist(lons):
    new_list = []
    oldval = 0
    # used to compare with the adjacent longitudes
    # and values exceed, disconnect linestring
    treshold = 10   
    for ix,ea in enumerate(lons):
        diff = oldval - ea
        if (ix>0):
            if (diff>treshold):
                ea = ea+360
        oldval = ea
        new_list.append(ea)
    return new_list

def mapfaults(ax,
              url='https://raw.githubusercontent.com/GEMScienceTools/gem-global-active-faults/master/geojson/gem_active_faults_harmonized.geojson',
              fallback_label='Active faults',
              linestyle_str = [(0, (3, 1, 1, 1, 1, 1)),':','-.','--'],
              legendfaults=True,
              **opt):
              
    if 'color' not in opt:
        opt['color']='r'
    if 'alpha' not in opt:
        opt['alpha']=.5
    if 'linewidth' not in opt:
        opt['linewidth']=0.8

    f=os.environ['HOME']+'/.local/share/cartopy/faults.geojson'
    if not glob.glob(f):
        print('downloading',url,'in',f)
        r = requests.get(url, allow_redirects=True)
        open(f, 'wb').write(r.content)

    faults=geopandas.read_file(f)
    extent=ax.get_extent(opt['transform'])
    labels=[]
    slip_types = list(set([slip_type for slip_type in faults.slip_type.values if None is not slip_type and  '_' in slip_type] ))

    for i,slip_type in enumerate(faults.slip_type.values):

        # grab x and y of the first geometry object
        geometry  = faults.iloc[i].geometry
        xs = geometry.xy[0].tolist()
        if not len( [x for x in xs if x>=extent[0] and x<=extent[1] ]):
            continue

        ys = geometry.xy[1].tolist()
        if not len( [x for x in ys if x>=extent[2] and x<=extent[3] ]):
            continue

        label=None

        # special cases
        tmpopt=opt.copy()
        if slip_type is not None and  '_' in slip_type:
            tmpopt['linewidth'] = opt['linewidth']*2
            tmpopt['linestyle']=linestyle_str[slip_types.index(slip_type)]
            if slip_type not in labels:
                label=slip_type.replace('_',' ')
                labels += [slip_type]

        elif fallback_label not in labels:
            label = fallback_label
            labels += [fallback_label]

        if not legendfaults:
            label=None

        # plot the geometry using sanitized values
        ax.plot(xs, ys, 
                label=label,
                path_effects=[matplotlib.patheffects.withStroke(linewidth=tmpopt['linewidth']*2,foreground="w",alpha=.5)], #
                **tmpopt)

class _TransformedBoundsLocator:
    """
    Axes locator for `.Axes.inset_axes` and similarly positioned Axes.
    The locator is a callable object used in `.Axes.set_aspect` to compute the
    axes location depending on the renderer.
    """

    def __init__(self, bounds, transform):
        """
        *bounds* (a ``[l, b, w, h]`` rectangle) and *transform* together
        specify the position of the inset Axes.
        """
        self._bounds = bounds
        self._transform = transform

    def __call__(self, ax, renderer):
        # Subtracting transSubfigure will typically rely on inverted(),
        # freezing the transform; thus, this needs to be delayed until draw
        # time as transSubfigure may otherwise change after this is evaluated.
        return mtransforms.TransformedBbox(
            mtransforms.Bbox.from_bounds(*self._bounds),
            self._transform - ax.figure.transSubfigure)

def makefigax(fig=None,ax=None,axprop={},figprop={}):
    if ax is not None:
        ax2 = ax.figure.add_axes(ax.get_position(True), 
                                **axprop, 
                                axes_locator=_TransformedBoundsLocator([0, 0, 1, 1], ax.transAxes))
        fig = ax2.figure
        ax.remove()
        ax=ax2
    else:
        if fig is not None:
            pass
        else :
            fig = matplotlib.pyplot.figure(**figprop)
        ax = matplotlib.pyplot.axes(**axprop)
    ax._autoscaleXon = False
    ax._autoscaleYon = False
    return ax

def bmap(bounds=[-89, -83, 8, 14],
         figsize=(5,5),
         rightside=True,
         top=False,
         fig=None,
         ax=None,
         legendfaults=True,
         bg=True,
         label_style={'size':'small',
                      'path_effects':[matplotlib.patheffects.withStroke(linewidth=2,foreground="w")]}):
    padlo=(bounds[1]-bounds[0])/10
    padla=(bounds[3]-bounds[2])/10
    midlo=(bounds[1]+bounds[0])/2
    mila=(bounds[3]+bounds[2])/2

    projection=ccrs.Orthographic(central_longitude=midlo,central_latitude=mila)
    
    ax=makefigax(fig=fig,
                 ax=ax,
                 axprop={'projection':projection},
                 figprop={'figsize':figsize})


    ax.set_extent([bounds[0]-padlo, bounds[1]+padlo,
                   bounds[2]-padla, bounds[3]+padla])
    
    gl=ax.gridlines(draw_labels=True, 
                    dms=True, 
                    x_inline=False, 
                    y_inline=False)
    
    if rightside:
        gl.right_labels= False
    else:
        gl.left_labels= False
    if top:
        gl.top_labels= False
    else:
        gl.bottom_labels= False

    gl.ylabel_style=label_style
    gl.xlabel_style=label_style
    gl.xpadding=7
    gl.ypadding=7

    ax.yaxis.tick_right()
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.xaxis.tick_top()
    ax.tick_params(axis="y", direction="out", length=4)
    ax.tick_params(axis="x", direction="out", length=4)


    #ax.add_feature(cfeature.LAND, alpha=0.95,color='w')
    #ax.add_feature(cfeature.OCEAN, alpha=0.95,color='w')
    #ax.add_feature(cfeature.LAKES, alpha=0.95,color='w')

    ax.add_feature(cfeature.BORDERS, linewidth=2,color='w')
    ax.add_feature(cfeature.BORDERS, linewidth=1,color='.5')
    ax.add_feature(cfeature.COASTLINE, linewidth=2,color='w')    
    ax.add_feature(cfeature.COASTLINE, linewidth=1,color='.5') 
    ax.add_feature(cfeature.RIVERS, linewidth=.5)
    
    if bg: 
        ax.add_image(cimgt.GoogleTiles(url='https://server.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}.jpg'), 8)

        ax.add_wms(wms='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv',
                layers=['GEBCO_LATEST'],
                alpha=1/3,
                zorder=1
                )
        if False:
            # Create a Stamen Terrain instance.
            stamen_terrain = cimgt.Stamen('terrain-background')
            # Add the Stamen data at zoom level 8.
            ax.add_image(stamen_terrain, 7,alpha=1/2)   
    
    
    mapcities(ax,transform=ccrs.Geodetic())
    mapfaults(ax,transform=ccrs.Geodetic(),legendfaults=legendfaults)
    mapcountrynames(ax,transform=ccrs.Geodetic())
    return ax

def _forward(x):
    return sqrt(x)
def _inverse(x):
    return x**2

def plot(Hyp,
         mtypes,
         ax,
         label='T$^{+}_{%s}$',
         times=None,
         labelmarkersize=64,
         facecolors={'Mfd':'C1','MVS':'C2','None':[1,0,1]},
         cmaps={'Mfd':'autumn','MVS':'winter'},
         nodt=False,
         **opt):
    
    norm = colors.FuncNorm((_forward, _inverse), 
                           vmin=0, 
                           vmax=40)

    scattermsizes=[]
    scattercmaps=[]

    if times is None:
        times=[xyz[4] for mtype in Hyp for xyz in Hyp[mtype]]    

    lima = [min(times),max(times)-min(times)]
    
    for mtype in Hyp:
        s=mtypes.index(mtype)
        x=[xyz[0]                            for xyz in Hyp[mtype]]
        y=[xyz[1]                            for xyz in Hyp[mtype]]
        z=[xyz[2]                            for xyz in Hyp[mtype]]
        m=[mag2size(xyz[3])                  for xyz in Hyp[mtype]]        
        # alpha by time
        a=[((xyz[4]-lima[0])/lima[1]+.5)*2/3 for xyz in Hyp[mtype]]       
        # alpha by depth
        #a=[((xyz[4]-lima[0])/lima[1]+.5)*2/3 for xyz in Hyp[mtype]]
        dt=[xyz[5]                           for xyz in Hyp[mtype]]

        if len(Hyp[mtype]) and len(Hyp[mtype][0])>6:
            xx=[[xyz[0],xyz[6]]              for xyz in Hyp[mtype]]
            yy=[[xyz[1],xyz[7]]              for xyz in Hyp[mtype]]
            dt=[xyz[11]                      for xyz in Hyp[mtype]]

        if nodt:
            dt=[ 0                           for xyz in Hyp[mtype]]

        o=argsort(m)
        
        mt=mtype[1:]
        if 'None' == mtype:
            mt=''
        addlab=' (%d ev.)'%len(m)

        ax.scatter([None], 
                   [None], 
                   [labelmarkersize], 
                   edgecolor='k', 
                   linewidths=.5,
                   facecolor=facecolors[mtype],#'C%d'%(len(mtypes)-1-s), 
                   label=label%mt+addlab, 
                   **opt)
        
        scattermsizes += [ax.scatter([None for i in o ], 
                              [None for i in o ], 
                              [m[i] for i in o], 
                              edgecolor='k', 
                              facecolor='0.5',
                              linewidths=.5,  
                              **opt)]
        if mtype in cmaps and mtype not in [im[1] for im in scattercmaps]:
            scattercmaps += [(ax.scatter([None], 
                                        [None], 
                                        [1],
                                        c=[1],
                                        cmap=cmaps[mtype],
                                        norm=norm),
                             "%s delay (s)"%mtype)]
        if False:
            scatters+=[ax.scatter([x[i] for i in o ], 
                                  [y[i] for i in o ], 
                                  [m[i] for i in o], 
                                  color='C%d'%(len(mtypes)-1-s), 
                                  linewidths=0.01, 
                                  label=label%mt+addlab, 
                                  **opt)]
            

        for j,i in enumerate(o):
            
            ax.scatter([x[i]], 
                       [y[i]], 
                       [m[i]], 
                       edgecolor='w', 
                       linewidths=4, 
                       alpha=a[i], 
                       **opt)
            
            ax.scatter([x[i]], 
                       [y[i]], 
                       [m[i]], 
                       edgecolor='k', 
                       linewidths=2, 
                       alpha=a[i], 
                       **opt)
            
            if mtype in cmaps:
                ax.scatter([x[i]], 
                        [y[i]], 
                        [m[i]],
                        c=[dt[i]],#color='C%d'%(len(mtypes)-1-s), 
                        cmap=cmaps[mtype],
                        norm=norm,
                        linewidths=0.01, 
                        alpha=a[i],
                        **opt)
            else:
                ax.scatter([x[i]], 
                        [y[i]], 
                        [m[i]],
                        color=facecolors[mtype],#'C%d'%(len(mtypes)-1-s), 
                        linewidths=0.01, 
                        alpha=a[i],
                        **opt)

            if len(Hyp[mtype]) and len(Hyp[mtype][0])>6:

                ax.plot(xx[i], 
                        yy[i],
                        alpha=a[i],
                        linewidth=2, 
                        color='w', 
                        solid_capstyle='round',
                        markersize=0, 
                        **opt)
                
                ax.plot(xx[i], 
                        yy[i],
                        alpha=a[i], 
                        linewidth=1, 
                        color=facecolors[mtype],#'C%d'%(len(mtypes)-1-s), #cmaps[mtypes](dt,norm=norm),#
                        solid_capstyle='round',
                        markersize=0,
                        **opt)

    return scattermsizes, scattercmaps
