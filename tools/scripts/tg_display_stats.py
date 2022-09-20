#!/usr/bin/python

#
# tg_display_stats
#
# ----------------------------------------------------------------------
# Copyright (c) 2011-2013 LEGOS/CTOH
# All rights reserved.
#
# Redistribution and use  in source  and binary  forms,  with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   1. Redistributions of  source  code must retain the  above copyright
#      notice, this list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright
#      notice,  this list of  conditions and the following disclaimer in
#      the  documentation  and/or  other   materials provided  with  the
#      distribution.
#
# THIS  SOFTWARE IS PROVIDED BY  THE  COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND  ANY  EXPRESS OR IMPLIED  WARRANTIES,  INCLUDING,  BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES  OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR  PURPOSE ARE DISCLAIMED. IN  NO EVENT SHALL THE COPYRIGHT
# HOLDERS OR      CONTRIBUTORS  BE LIABLE FOR   ANY    DIRECT, INDIRECT,
# INCIDENTAL,  SPECIAL,  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN  CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE   OF THIS SOFTWARE, EVEN   IF ADVISED OF   THE POSSIBILITY OF SUCH
# DAMAGE.
#
# DATE OF CREATION: 2012-08-06
# AUTHOR: Sara Fleury  sara.fleury@legos.obs-mip.fr
# CONTACT: ctoh_products@legos.obs-mip.fr
# ----------------------------------------------------------------------


"""
DESCRIPTION: Display stats on tides at stations

USAGE:  tg_display_stats.py <file*.stat>  [OPTIONS]
or   :  tg_display_stats.py "file*"       [OPTIONS] 

OPTIONS:
      -l                : print list of station
      -m[<quality>]     : plot stations on a map, with shoreline_quality in %s
      -p[abn]           : plot data a (after detide) ; and/or b (before detide) ; and/or n (no label)
      -v<level>         : verbose level

EXAMPLES:
     ./tg_display_stats.py /gsa/fles/gloss/hourly-20130416/analysis.ondes65.FES2004.01y/h82*.stat -ml
    ./tg_display_stats.py /gsa/fles/gloss/hourly-20130416/analysis.ondes65.FES2004.01y/*.stat -l
    ./tg_display_stats.py /gsa/fles/gloss/hourly-20130416/analysis.ondes65.FES2004.01y/h826_STOCKHOL*.stat -v3
    ./tg_display_stats.py /gsa/fles/gloss/hourly-20130416/analysis.ondes65.FES2004.01y/h8*.stat -pab
    ./tg_display_stats.py /gsa/fles/gloss/hourly-20130416/analysis.ondes65.FES2004.01y/h8*.stat -m -ql

"""


#------------------------------------------------------------------#
# External modules
#------------------------------------------------------------------#

# generic
import os,sys,math
import glob
import numpy as np
import inspect



# ----------------------------------------------------------------------
# usage
#

def stat_error(str):
    #sys.stdout = sys.stderr
    print ("** ERROR %s:%s: %s" % (__file__,inspect.currentframe().f_back.f_lineno, str))
    print ""
    sys.exit()

def stat_warning(str):
    #sys.stdout = sys.stderr
    print ("!! WARNING %s:%s: %s" % (__file__,inspect.currentframe().f_back.f_lineno, str))
    sys.stdout = sys.__stdout__

def stat_info(str):
    #sys.stdout = sys.stderr
    print ("%s: %s" % (os.path.basename(__file__), str))

def stat_usage(str):
    #sys.stdout = sys.stderr
    print __doc__ % (BASEMAP_QUALITY)
    print ("** ERROR %s:%s: %s\n" % (__file__,inspect.currentframe().f_back.f_lineno, str))
    sys.exit()



#------------------------------------------------------------------#
# plot modules
#------------------------------------------------------------------#
global PLOT_FLAG
PLOT_FLAG=None
try:
    import matplotlib.pyplot as plt
    PLOT_FLAG="1d"
    try:
        from mpl_toolkits.basemap import Basemap
        PLOT_FLAG="2d"
    except:
        stat_warning ("cannot plot map without module 'matplotlib'")

except:
    stat_warning ("cannot plot without python module 'matplotlib'")
    PLOT_FLAG = None


# ----------------------------------------------------------------------
# plot coastlines and so on
global BASEMAP_RESOLUTION
def import_basemap ():
    global BASEMAP_RESOLUTION
    if BASEMAP_RESOLUTION:
        return
    BASEMAP_RESOLUTION  = {'c':10000,'l':1000,'i':100,'h':10,'f':1}

BASEMAP_QUALITY = {'c':'coarse (res=10000)','l':'low (res=1000)','i':'intermediate (res=100)','h':'hight (res=10)','f':'full (res=1)'}


#------------------------------------------------------------------#
# Read stat files
#------------------------------------------------------------------#

def read_stat (f_list_in, stat_flag, verbose=1):

    stations_stat = {}
    stations_list = {}

    if len(f_list_in) == 1 and '*' in f_list_in[0]:
        f_list_in += '.stat'
        f_list = glob.glob(f_list_in[0])
    else:
        f_list = f_list_in

    if not f_list:
        stat_error("No files %s" % f_list_in)
    f_list.sort()

    for f_full in f_list:

        f = os.path.basename(f_full)

        # decode name, dates
        full_name = f[:14]
        code = f[1:4]
        name = f[5:14]
        station = f[:4]
        pos_tides = f[-9].rfind('.')
        tides = f[pos_tides:-9]
        waves = f[f[pos_tides].rfind('.'):pos_tides]
        start = f[14:22]
        end = f[23:31]

        # record station
        if station in stations_list.keys() and  not stat_flag:
            continue
            
        # open input ascii file
        if verbose > 3:
            stat_info( "read_stat: Read %s" % f)
        f_in = open(f_full)
        lines = f_in.readlines()
        f_in.close()
        # (nb_samp, mean, var, dev) = lines[2]
        pos = (lines[0].split())[-2:]
        pos[0] = float(pos[0])
        pos[1] = float(pos[1])

        if not station in stations_list.keys():
            stations_list[station] = (name, code, pos)
            if not stat_flag:
                continue
            
        stat_post_detide = lines[3].split()
        stat_pre_detide = lines[2].split()
        var_dev = (float(stat_pre_detide[3]) - float(stat_post_detide[3]))/float(stat_pre_detide[3]) *100.

        # record station
        if not full_name in stations_stat:
            stations_stat[full_name] = (name, code, pos, [])
        stations_stat[full_name][3].append((start, end, stat_pre_detide, stat_post_detide, var_dev))

    # end
    #stations_list.sort()
    return stations_list, stations_stat

#------------------------------------------------------------------------------#
# List stats
#------------------------------------------------------------------------------#

def list_stations (stations_list, verbose=1): 

    print "%d  XYN" % len(stations_list)
    sorted_list = stations_list.keys()
    sorted_list.sort()
    for s in sorted_list:
        v = stations_list[s]
        print "%9.3f %9.3f   %s_%s" % (v[2][0],v[2][1],s,v[0])
        #print "%s_%s" % (s,v[0])

#------------------------------------------------------------------------------#
# Print stats
#------------------------------------------------------------------------------#

def print_stat (stations_stat, verbose=1):

    head ="""
# -----------------------------------------------------------------------------
# CODE_STATION         POS             NB       DEV_REDUCTION         DEV (m) 
#                  lon      lat       years   worst  best  avrg      worst avrg
# -----------------------------------------------------------------------------
%d  NXY"""

    
    stations_list = stations_stat.keys()
    stations_list.sort()
    print head % len(stations_list)


    if verbose == 1:
        print "#  start   dev  var_dev    (after de-tide)"
    if verbose == 2:
        print "#  start   end   :  nb_samp  mean    var     dev   var_dev   (after de-tide)"
    if verbose == 3:
        print "#  start   end   :  nb_samp  mean    var     dev             (before de-tide)"
        print "#                   nb_samp  mean    var     dev   var_dev   (after  de-tide)"

    for s in stations_list:
        nb_files=0
        devs = []
        var_devs = []
        v=stations_stat[s]
        if verbose > 0:
            print "%s %8.3f  %8.3f  " % (s,v[2][0],v[2][1])
        periods = v[3]
        for (start, end, stat_pre_detide, stat_post_detide, var_dev) in periods:
            nb_files += 1
            var_devs.append(var_dev)
            devs.append(float(stat_post_detide[-1]))
            
            if verbose == 1:
                (nb_samp, mean, var, dev) = stat_post_detide
                print "%8s %6s %2.0f%%" % (start, dev, var_dev)
                
            elif verbose == 2:
                (nb_samp, mean, var, dev) = stat_post_detide
                print "%8s %8s : %5s %8s  %6s  %6s %2.0f%%" % (start, end,  nb_samp, mean, var, dev, var_dev)
                
            elif verbose == 3:
                (nb_samp, mean, var, dev) = stat_pre_detide
                print "%8s %8s : %5s %8s  %6s  %6s" % (start, end,  nb_samp, mean, var, dev)
                (nb_samp, mean, var, dev) = stat_post_detide
                print "%8s %8s : %5s %8s  %6s  %6s %2.0f%%" % ("", "", nb_samp, mean, var, dev, var_dev)

        if verbose >= 0:
            print "%s   %8.3f  %8.3f  " % (s,v[2][0],v[2][1]),
        else:
            print "%s" % (s),

        # print synthesis per station
        best_var = max(var_devs)
        worst_var = min(var_devs)
        average_var = np.average(var_devs)
        worst_dev = max(devs)
        average_dev = np.average(devs)
        if np.isnan(best_var): best_var_str = " nan"
        else: best_var_str = "%3d%%" % int(best_var)
        if np.isnan(worst_var): worst_var_str = " nan"
        else: worst_var_str = "%3d%%" % int(worst_var)
        if np.isnan(average_var): average_var_str = " nan"
        else: average_var_str = "% 3d%%" % int(average_var)
        if np.isnan(worst_dev): worst_dev_str = " nan"
        else: worst_dev_str = "%3d%%" % int(worst_dev)
        if np.isnan(average_dev): average_dev_str = " nan"
        else: average_dev_str = "% 3d%%" % int(average_dev)
        
        print "%3d   %3s  %3s   %3s   %7.3f %7.3f" % (nb_files, worst_var_str, best_var_str, average_var_str, worst_dev, average_dev)

        
#------------------------------------------------------------------------------#
# Print stats
#------------------------------------------------------------------------------#

def plot_stat (stations_stat, plot_type='a', verbose=1):

    global PLOT_FLAG
    if not PLOT_FLAG:
        stat_warning("Cannot plot on the machin")
        return

    COLOR_TAB = ['b','r','g','c','m','y','k',
                 'orange', 'grey',   
                 'lightgreen', 'pink',
                 'teal', 'darkred', 'darkgreen', 'darkblue', 'darkorange']

    if plot_type == False:
        plot_type = 'a'
    elif not a in plot_type and not 'b' in plot_type:
        plot_type += 'a'
        
    stations_list = stations_stat.keys()
    stations_list.sort()
    y_plot_style = '.'
    col = 0

    fig = plt.figure()
    ax = fig.gca()
    legend = ''

   
    for s in stations_list:
        devs_pre = []
        devs_post = []
        time = []
        
        v=stations_stat[s]
        slabel = "%15s %6.3f %6.3f" % (s, v[2][0], v[2][1])
        periods = v[3]
        for (start, end, stat_pre_detide, stat_post_detide, var_dev) in periods:
            devs_pre.append(float(stat_pre_detide[-1]))
            devs_post.append(float(stat_post_detide[-1]))
            time.append(int(start[:4]))

        if 'b' in plot_type:
            if not 'a' in plot_type and not 'n' in plot_type:
                plt.plot(time, devs_pre, '+', color=COLOR_TAB[col], label=slabel)
            else:
                plt.plot(time, devs_pre, '+', color=COLOR_TAB[col])
            plt.plot(time, devs_pre, '--', color=COLOR_TAB[col])
        if 'a' in plot_type:
            if not 'n' in plot_type:
                plt.plot(time, devs_post, '.', color=COLOR_TAB[col], label=slabel)
            else:
                plt.plot(time, devs_post, '.', color=COLOR_TAB[col])
            plt.plot(time, devs_post, '-', color=COLOR_TAB[col])
        col += 1
        if col >= len(COLOR_TAB): col=0
        
    plt.ylabel("dev (m)", fontsize='large')
    plt.xlabel("year")
    title = "Deviation "
    if 'b' in plot_type: title += "before (- -) "
    if 'a' in plot_type: title += "after (___) "
    title += "detiding"
    
    plt_title = plt.title(title)
    plt_title.set_fontsize('large')
    plt.grid(True)
    plt.minorticks_on()

    if not 'n' in plot_type:
        leg = plt.legend(labelspacing=0,loc='best')
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
        
    # plot
    try:
        plt.show()
    except:
        sys.exit()
            
#------------------------------------------------------------------#
# Plot map
#------------------------------------------------------------------#

def plot_map (stations_list, shoreline_quality='l', verbose=1):

    global PLOT_FLAG
    if PLOT_FLAG != "2d":
        stat_warning("Cannot plot map on the machin")
        return

    import_basemap()
    
    lons = []
    lats= []
    nb_samps = []
    nb_ys = []
    var_devs = []
    file_name = None
    
    # read stats
    for s,v in stations_list.iteritems():
        lon = v[2][0]
        lat = v[2][1]
        lons.append(lon)
        lats.append(lat)

    # map size
    lonmin = min(lons)
    lonmax = max(lons)
    latmin = min(lats)
    latmax = max(lats)
    lonmin -= 10
    lonmax += 10
    latmin -= 10
    latmax += 10
    if latmin < -90: latmin = -90
    if latmax > 90: latmax = 90
    
    
    # create figure
    fig = plt.figure(1) # , facecolor = 'w', bgcolor='w')
    mapframe = plot_shoreline(fig, lonmin, lonmax, latmin, latmax, shoreline_quality)

    # add stations to shorelines
    x, y = mapframe(lons,lats)
    plt.scatter(x,y,marker='+', s=20, color='r')

    for s,v in stations_list.iteritems():
        name  = s + '_' + v[0]
        lon = v[2][0]
        lat = v[2][1]
        plt.annotate(
            name, 
            xy = (mapframe(lon,lat)), xytext = (-0, 0),
            size='small',
            color='r',
            textcoords = 'offset points', ha = 'right', va = 'bottom')
        
    plt.title('Locations of %s GLOSS stations' % len(lons))

    # show all

    # on file
    #file_name='stations.png'
    if file_name:
        plt.savefig(file_name)
        stat_info( "display %s" % file_name)
        sys.exit()
    # on screen
    try:
        plt.show()
    except:
        sys.exit()

  
#------------------------------------------------------
# 
#
# where to put meridians and parallels labels
label_merid = [1, 0, 0, 1]
label_paral = [1, 0, 0, 1]

def plot_shoreline (fig, lon_min, lon_max, lat_min, lat_max, shoreline_quality):

    # minimum space between meridians and parallels: 10 deg
    grid_res_deg = [10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    grid_nb_min = 4

    # minimum number of drawn meridians and parallels
    dX = lon_max - lon_min
    if dX < 0: dX=dX+360
    for res in grid_res_deg:
        if dX/res >= grid_nb_min:
            break
    d_merid = res


    dY = lat_max - lat_min
    for res in grid_res_deg:
        if dY/res >= grid_nb_min:
            break
    d_paral = res
    
    parallels = np.arange(-90, 90.+d_paral, d_paral).tolist()
    meridians = np.arange(-180, 360, d_merid).tolist()

    # shoreline
    if not shoreline_quality:
        shoreline_quality = 'l'

    # main frame
    #if lon_min > 180 or lon_max > 180:
    #    lon_min = lon_min - 360
    #    lon_max = lon_max - 360
    mapframe = Basemap(llcrnrlon = lon_min, llcrnrlat = lat_min, urcrnrlon = lon_max, urcrnrlat = lat_max, resolution = shoreline_quality) # , area_thresh = 1.)

    # set in grey _FillValue
    #palette  =  plt.cm.jet
    #palette.set_bad('gray',  1.0)

    # parallels and meridians
    mapframe.drawparallels(parallels, labels = label_paral, fontsize = 'small')
    mapframe.drawmeridians(meridians, labels = label_merid, fontsize = 'small')

    # coastlines
    mapframe.drawcoastlines()
    #mapframe.drawrivers(linewidth=0.01, color="lightgrey", #030303
    #                    antialiased=1, ax=None, zorder=1)
    #mapframe.fillcontinents()

    # display data color map
    #data = np.where(np.isnan(data), 1.e10, data)
    #data[np.abs(data)>900] = 1.e10
    #data =  np.ma.masked_values(data,  1.e10)

    # show with no interpolation
    ##im = mapframe.imshow(np.squeeze(data), palette, norm, interpolation = 'nearest') #, interpolation = 'hanning')
    #fig = plt.gcf()
    #ax = plt.gca()

    # draw
    #fig.canvas.draw()

    # show all
    #plt.show()

    return mapframe


 
# ----------------------------------------------------------------------
# main
#

if __name__ == '__main__':
    global BASEMAP_RESOLUTION
    BASEMAP_RESOLUTION = None

    verbose = 0
    shoreline_quality = 'i'

    # default options
    stat_flag = True
    list_flag = False
    plot_flag = False
    map_flag = False
    plot_type = "a"
    f_list = []
    
    if len(sys.argv) < 2:
        stat_usage("Bad params")

    for a in sys.argv[1:]:
        if a[0] == '-':
            for o in a[1:]:
                if o == 'v':
                    verbose = int(a[2])
                elif o == 'l':
                    list_flag = True
                    stat_flag = False
                elif o[0] == 'm':
                    map_flag = True
                    list_flag = True
                    stat_flag = False
                    if len(o) == 2:
                        shoreline_quality = o[1]
                elif o == 'p':
                    plot_flag = True
                    if len(a)>2:
                        plot_type = a[2:]
                        for c in plot_type:
                            if c != 'a' and c != 'b' and c != 'n':
                                stat_usage("Plot options are: a (after detide) ; and/or b (before detide) ; and/or n (no label)")

                #else:
                #    stat_usage("Unknown option %s" % a)
        else:
            f_list.append(a)
            
    # read files
    stations_list, stations_stat = read_stat(f_list, stat_flag, verbose)

    # list stations
    if list_flag:
        list_stations(stations_list, verbose)

    # stats
    if stat_flag:
        print_stat(stations_stat, verbose)
        
    # plot map
    if map_flag:
        plot_map(stations_list, shoreline_quality, verbose)

    # plot stat
    if plot_flag:
        plot_stat(stations_stat, plot_type, verbose)
        
   

# WORST STATIONS (dev after detiding > 0.3)
#
# h906_MOULMEIN_     97.617    16.483      6    78%   53%    62%     0.551   0.430
# h825_CUXHAVEN_      8.717    53.867     24    69%   54%    64%     0.499   0.397
# h595_NOME_____   -165.430    64.500     19     6%    3%     4%     0.427   0.364
# h124_CHITTAGO_     91.827    22.237      6    85%   75%    79%     0.334   0.283
# h336_BOOBY_IS_    141.917   -10.600     25    80%   60%    69%     0.331   0.259
# h601_ESPERANZ_    -56.995   -63.395     13    86%   42%    77%     0.396   0.152
# h729_MAR_DEL__    -57.538   -38.035      9    48%   31%    36%     0.300   0.269

# 2nd WORST STATIONS (dev after detiding > 0.2)
#
# h007_MALAKAL__    134.463     7.330     33    91%   58%    78%     0.204   0.100
# h049_MINAMITO_    153.967    24.300     16    13%    6%    10%     0.288   0.237
# h125_PRIGI____    111.733    -8.283      6    74%   61%    67%     0.217   0.179
# h128_THEVENAR_    133.640   -32.150     15    62%   49%    54%     0.235   0.209
# h155_DZAOUDZI_     45.257   -12.780      5    95%   71%    87%     0.235   0.109
# h162_CILACAP__    109.017    -7.752      6    69%   56%    63%     0.203   0.166
# h186_KNYSNA___     23.050   -34.083     20    81%   58%    69%     0.209   0.142
# h188_RICHARD'_     32.083   -28.800     24    95%   50%    78%     0.255   0.109
# h210_FLORES___    -31.120    39.453     21    80%   40%    66%     0.218   0.107
# h260_DUCK_PIE_    -75.740    36.183     33    63%   50%    58%     0.203   0.167
# h261_CHARLEST_    -79.933    32.783     33    76%   64%    72%     0.212   0.168
# h264_ATLANTIC_    -74.418    39.355     33    67%   57%    63%     0.203   0.177
# h273_BASQUES__    -59.133    47.567     16    67%   37%    62%     0.232   0.138
# h274_CHURCHIL_    -94.200    58.783     33    85%   75%    81%     0.281   0.216
# h281_CANANEIA_    -47.925   -25.017     28    57%   44%    48%     0.226   0.205
# h295_STORNOWA_     -6.388    58.208     31    89%   80%    83%     0.217   0.182
# h302_BALBOA___    -79.573     8.962     31    91%   84%    89%     0.220   0.151
# h328_KO_LAK___     99.817    11.795     26    63%   55%    60%     0.212   0.198
# h383_VUNG_TAU_    107.072    10.340      6    89%   75%    81%     0.205   0.155
# h399_LORD-HOW_    159.067   -31.533     16    76%   59%    66%     0.205   0.168
# h417_SADENG___    110.783    -8.500      5    70%   58%    65%     0.220   0.178
# h540_PRINCE_R_   -130.333    54.317     33    90%   85%    88%     0.226   0.183
# h542_TOFINO___   -125.917    49.150     33    86%   73%    80%     0.226   0.166
# h558_NEAH_BAY_   -124.617    48.368     33    84%   70%    77%     0.220   0.168
# h559_SITKA____   -135.342    57.052     33    85%   77%    82%     0.216   0.167
# h560_SEWARD___   -149.427    60.120     33    84%   76%    81%     0.243   0.187
# h570_YAKUTAT__   -139.735    59.547     33    84%   77%    81%     0.216   0.180
# h571_KETCHIKA_   -131.625    55.333     33    90%   85%    87%     0.218   0.182
# h574_SAND_POI_   -160.502    55.337     17    77%   65%    72%     0.250   0.196
# h579_PRUDHOE__   -148.527    70.400     20     8%    2%     4%     0.259   0.186
# h592_SOUTH_BE_   -124.043    44.625     33    85%   72%    79%     0.216   0.159
# h599_DIEGO_RA_    -68.715   -56.508      8    60%   50%    56%     0.200   0.167
# h600_USHUAIA__    -68.295   -54.805      4    73%   58%    63%     0.204   0.177
# h684_PUERTO_M_    -72.967   -41.483     33    93%   86%    90%     0.210   0.147
# h731_PUERTO_M_    -65.030   -42.763      3    84%   80%    82%     0.283   0.254
# h752_FORT_PUL_    -80.902    32.033     33    79%   70%    76%     0.227   0.183
# h767_GALVESTO_    -94.788    29.285     32    38%   21%    32%     0.216   0.179
# h775_GALVESTO_    -94.793    29.310     33    28%   14%    23%     0.209   0.172
# h800_ANDENES__     16.150    69.317     22    74%   62%    67%     0.208   0.177
# h803_RORVIK___     11.250    64.867     33    76%   65%    70%     0.218   0.185
# h805_VARDO____     31.100    70.333     29    81%   71%    77%     0.229   0.175
# h809_SCORESBY_    -21.983    70.483      6    57%   35%    45%     0.204   0.158
# h818_SMOGEN___     11.215    58.350     13    11%    7%     9%     0.232   0.196
# h819_GOTHENBU_     11.800    57.683     33     5%    2%     3%     0.247   0.202
# h826_STOCKHOL_     18.083    59.317     32     0%    0%     0%     0.260   0.193
# h833_NAIN_____    -61.700    56.550     12    80%   62%    74%     0.215   0.154
# h834_MALIN_HE_     -7.333    55.367      5    81%   73%    78%     0.220   0.184
# h907_SITTWE(A_     92.900    20.150      4    73%   62%    69%     0.246   0.198
# h908_PORT_BLA_     92.767    11.683      2    82%   56%    69%     0.223   0.155
