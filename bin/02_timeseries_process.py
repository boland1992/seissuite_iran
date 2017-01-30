#!/usr/bin/python -u
"""
[Advice: run this script using python with unbuffered output:
`python -u 01_timeseries_process.py`]

This script reads seismic waveform data from a set of stations, and
calculates the cross-correlations between all pairs of stations. The
data (in miniseed format) must be located in folder *MSEED_DIR*. The
stations information (coordinates, instrument response) can be read
from dataless seed files (if *USE_DATALESSPAZ* = True) located in
folder *DATALESS_DIR*, and/or stationXML files (if *USE_STATIONXML* =
True) located in folder *STATIONXML_DIR*. Note that two different
stations MUST HAVE DIFFERENT NAMES, even if they do not belong to
the same network. Also, one given station cannot have several
sets of coordinates: if so, it will be skipped.

In the current version of the program, miniseed files MUST be
organized inside their directory as:
<year>-<month>/<network>.<station>.<channel>.mseed, e.g.:
1988-10/BL.JFOB.BHZ.mseed
So, there is one sub-directory per month, and inside it, one miniseed
file  per month and per station.

The implemented algorithm follows the lines of Bensen et al.,
"Processing seismic ambient noise data to obtain reliable broad-band
surface wave dispersion measurements", Geophys. J. Int. (2007).

The procedure consists in stacking daily cross-correlations between
pairs of stations, from *FIRSTDAY* to *LASTDAY* and, in each given day,
rejecting stations whose data fill is < *MINFILL*. Define a subset of
stations to cross-correlate in *CROSSCORR_STATIONS_SUBSET* (or let it
empty to cross-correlate all stations). Define a list of locations to
skip in *CROSSCORR_SKIPLOCS*, if any. The cross-correlations are
calculated between -/+ *CROSSCORR_TMAX* seconds.

Several pre-processing steps are applied to the daily seismic waveform
data, before the daily cross-correlation is calculated and stacked:

(1) removal of the instrument response, the mean and the trend;

(2) band-pass filter between *PERIODMIN* and *PERIODMAX* sec

(3) down-sampling to sampling step = *PERIOD_RESAMPLE* sec

(4) time-normalization:

    - if *ONEBIT_NORM* = False, normalization of the signal by its
      (smoothed) absolute amplitude in the earthquake period band,
      defined as *PERIODMIN_EARTHQUAKE* - *PERIODMIN_EARTHQUAKE* sec.
      The smoothing window is *PERIODMAX_EARTHQUAKE* / 2;

    - if *ONEBIT_NORM* = False, one-bit normalization, wherein
      only the sign of the signal is kept (+1 or -1);

(5) spectral whitening of the Fourier amplitude spectrum: the Fourier
    amplitude spectrum of the signal is divided by a smoothed version
    of itself. The smoonthing window is *WINDOW_FREQ*.

Note that all the parameters mentioned above are defined in the
configuration file.

When all the cross-correlations are calculated, the script exports
several files in dir *CROSS*

"""
import os
import warnings
import datetime as dt
import itertools as it
import pickle
import obspy.signal.cross_correlation
import time
import glob
import sqlite3 as lite
import shutil
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import Trace


# set epoch timestamp 
epoch = dt.datetime(1970, 1, 1)

total_verbose = True
psd = False


# DECLUSTER STATIONS!
# remove stations that are too close to one another (set by degree radius!)
#from seissuite.spacing.search_station import Coordinates

#import matplotlib.pyplot as plt
#import numpy as np
# turn on multiprocessing to get one merged trace per station?
# to preprocess trace? to stack cross-correlations?
MULTIPROCESSING = {'merge trace': True,
                   'process trace': True,
                   'cross-corr': True}
# how many concurrent processes? (set None to let multiprocessing module decide)
NB_PROCESSES = None
if any(MULTIPROCESSING.values()):
    import multiprocessing as mp

# create a list of configuration files that will be iterated over! 
# must be 1 or more
from seissuite.ant.psconfig import (create_config_list, run_config, 
                                    remove_config)

config_list = create_config_list()
config_file = config_list[0]

total_time0 = dt.datetime.now()
#for config_file in config_list:
for i in range(0, 100):   
    
    # global variables MUST be defined 
    # with the function in the seissuite.ant.psconfig module 
    run_config(config_file)
  
  
    from seissuite.ant import (pscrosscorr, psstation, pspreprocess, pserrors, 
                               psstationSQL)

    # import CONFIG class initalised in ./configs/tmp_config.pickle
    config_pickle = 'configs/tmp_config.pickle'
    f = open(name=config_pickle, mode='rb')
    CONFIG = pickle.load(f)
    f.close()
    
    # import variables from initialised CONFIG class.
    MSEED_DIR = CONFIG.MSEED_DIR
    DATABASE_DIR = CONFIG.DATABASE_DIR
    DATALESS_DIR = CONFIG.DATALESS_DIR
    STATIONXML_DIR = CONFIG.STATIONXML_DIR
    CROSSCORR_DIR = CONFIG.CROSSCORR_DIR
    USE_DATALESSPAZ = CONFIG.USE_DATALESSPAZ
    USE_STATIONXML = CONFIG.USE_STATIONXML
    CROSSCORR_STATIONS_SUBSET = CONFIG.CROSSCORR_STATIONS_SUBSET
    CROSSCORR_SKIPLOCS = CONFIG.CROSSCORR_SKIPLOCS
    FIRSTDAY = CONFIG.FIRSTDAY
    LASTDAY = CONFIG.LASTDAY
    MINFILL = CONFIG.MINFILL
    FREQMIN = CONFIG.FREQMIN
    FREQMAX = CONFIG.FREQMAX
    CORNERS = CONFIG.CORNERS
    ZEROPHASE = CONFIG.ZEROPHASE
    PERIOD_RESAMPLE = CONFIG.PERIOD_RESAMPLE
    ONEBIT_NORM = CONFIG.ONEBIT_NORM
    FREQMIN_EARTHQUAKE = CONFIG.FREQMIN_EARTHQUAKE
    FREQMAX_EARTHQUAKE = CONFIG.FREQMAX_EARTHQUAKE
    WINDOW_TIME = CONFIG.WINDOW_TIME
    WINDOW_FREQ = CONFIG.WINDOW_FREQ
    XCORR_INTERVAL = CONFIG.XCORR_INTERVAL
    CROSSCORR_TMAX = CONFIG.CROSSCORR_TMAX
    PLOT_CLASSIC = CONFIG.PLOT_CLASSIC
    PLOT_DISTANCE = CONFIG.PLOT_DISTANCE
    MAX_DISTANCE = CONFIG.MAX_DISTANCE
    RANDOM_STACK = False
    RESP_REMOVE = False
    
    #FULL_COMB = CONFIG.FULL_COMB

    # initialise the required databases if they haven't already been.
    #if no two SQL databases exist, then create them! 
    TIMELINE_DB = os.path.join(DATABASE_DIR, 'timeline.db')
    RESP_DB = os.path.join(DATABASE_DIR, 'response.db')

   # RESP_DB = os.path.join(DATABASE_DIR, 'response.db')
    
   # if not os.path.exists(RESP_DB):
        # initialise response database for use with automated data selection!

  #      lite.connect(RESP_DB)
  #      from seissuite.database import response_database
    print TIMELINE_DB
    if not os.path.exists(RESP_DB):
        try:
            lite.connect(RESP_DB)
            print "\nCreating response database. Please be patient ... "
            from seissuite.database import response_database
        except:
            print "Response database could not be initialised ... "

    if not os.path.exists(TIMELINE_DB):
        # initialise timeline database to help the application find files!
        lite.connect(TIMELINE_DB)
        print "\nCreating timeline database. Please be patient ... "
        from seissuite.database import create_database                                        
    
    
    if psd:
        import powerdensity
        

    print "\nProcessing parameters:"
    print "- dir of miniseed data: " + MSEED_DIR
    print "- dir of dataless seed data: " + DATALESS_DIR
    print "- dir of stationXML data: " + STATIONXML_DIR
    print "- output dir: " + CROSSCORR_DIR
    print "- cross-correlation length (mins): " + str(XCORR_INTERVAL)
    print "- cross-correlation maximum time interval (s): " + str(CROSSCORR_TMAX)
    
    
    print "- band-pass: {:.1f}-{:.1f} s".format(1.0 / FREQMAX, 1.0 / FREQMIN)
    if ONEBIT_NORM:
        print "- normalization in time-domain: one-bit normalization"
    else:
        s = ("- normalisation in time-domain: "
             "running normalisation in earthquake band ({:.1f}-{:.1f} s)")
        print s.format(1.0 / FREQMAX_EARTHQUAKE, 1.0 / FREQMIN_EARTHQUAKE)
    fmt = '%d/%m/%Y'
    s = "- cross-correlation will be stacked between {}-{}"
    print s.format(FIRSTDAY.strftime(fmt), LASTDAY.strftime(fmt))
    subset = CROSSCORR_STATIONS_SUBSET
    if subset:
        print "  for stations: {}".format(', '.join(subset))
    print
    
        
    # Initializing collection of cross-correlations
    xc = pscrosscorr.CrossCorrelationCollection()
    
    #create a metadata list, may need dictionary based on how much info required         
    metadata = [] 
    
    #ask if system has crashed or stopped before another process was finished?
    
    print "\nScanning for partial pickle cross-correlation files ... "
    
    #maybe create pause statement for interesting load menu.
    
    # loading cross-correlations (looking for *.part.pickle files in folders in
    #in dir *CROSSCORR_DIR*)
    folder_list = sorted(glob.glob(os.path.join(CROSSCORR_DIR, '*')))
    
    pickle_list = []
    index = 0 #creating index for call 
    for folder in folder_list:
        #check to see if there are any pickle files in the xcorr time folder 
        if len(glob.glob(os.path.join(folder, '*.part.pickle'))) < 1:
            #print("There are no .pickle files in this folder. Skipping ...")
            continue
        else:
            #append name of pickle file path location string to pickle_list 
            pickle_list.append(glob.glob(os.path.join(folder, \
            '*.part.pickle'))[0])
    
    if len(pickle_list) < 1: 
        print("\nThere are no partial pickle files to begin again from.")
        print("\nThe program will start from the beginning")
        res = ""
        
    else:
        print "\nPlease choose a file to begin again from, or a combination thereof."
        print "Else hit enter to continue anew"
        #print combinations of partial pickle files available
        print '\n0 - All except backups (*~)'    
        print '\n'.join('{} - {}'.format(i + 1, f.split('/')[-2])
            for i, f in enumerate(pickle_list))
                
                
        #change folder_list to pickle_list if this gives problems
        #res = False#raw_input('\n')
        res = ''
        
        
        #res = raw_input('\n')



    
    #IF LIST INDEX OUT OF RANGE START PROGRAM ALSO    
    
    #if beginning again, reset time-series intervals to the where the selected 
    # .part.pickle file left off! 
    
    if not res:
        # ========================================
        #set output file name as normal
        # ========================================
        time_string = str(time.strftime("%d.%m.%Y") + "-" + time.strftime("%X"))
        responsefrom = []
        if USE_DATALESSPAZ:
            responsefrom.append('datalesspaz')
        if USE_STATIONXML:
            responsefrom.append('xmlresponse')
        OUTBASENAME_PARTS = [
        'XCORR-STACK',
        '-'.join(s for s in CROSSCORR_STATIONS_SUBSET) \
        if CROSSCORR_STATIONS_SUBSET else None,
        '{}-{}'.format(FIRSTDAY.strftime("%d.%m.%Y"),
                       LASTDAY.strftime("%d.%m.%Y")),
        '1bitnorm' if ONEBIT_NORM else None,
        '+'.join(responsefrom)
        ]
        
        OUTFILESNAME = '_'.join(p for p in OUTBASENAME_PARTS if p)
    
        OUTFILESPATH = os.path.join(CROSSCORR_DIR, time_string, OUTFILESNAME)
    
        OUTFOLDERS = os.path.join(CROSSCORR_DIR, 
                                  time_string,  
                                  'XCORR_PLOTS')

        
        OUT_SNR = os.path.join(CROSSCORR_DIR, time_string,  'SNR_PLOTS')
                                  
        #create unique folder in CROSS output folder named by the present time.
        if not os.path.exists(OUTFOLDERS):\
        os.makedirs(OUTFOLDERS)
        
        if not os.path.exists(OUT_SNR):\
        os.makedirs(OUT_SNR)
        
        # copy configuration file to output so parameters are known for each run                         
        OUTCONFIG = os.path.join(CROSSCORR_DIR, time_string, 
                                 os.path.basename(config_file))
        
        print 'Copying configuration file to output directory ... ' 
        shutil.copy(config_file, OUTCONFIG)    
        

        METADATA_PATH = '{}metadata.pickle'.format(OUTFILESPATH.\
                  replace(os.path.basename(OUTFILESPATH), ""))
    
    else:
        # ========================================
        #reset time as previous time, reset output paths as previous path name
        #reset cross-correlation dictionaries 
        # ========================================
        print 
        PART_PICKLE = pickle_list[int(res)-1]
        OUTFILESPATH = PART_PICKLE[:-12]
        out_basename = os.path.basename(OUTFILESPATH)
        
        
        print "Opening {} partial file for restart ... ".format(out_basename)

        # re-initialising .part.pickle collection of cross-correlations
        xc = pscrosscorr.load_pickled_xcorr(PART_PICKLE)
        
                    
        for key in xc.keys():
            for key2 in xc[key].keys():
                #help(xc[key][key2])
                #print xc[key][key2].endday
                a=5
                
                
        #most recent day last endday of list
        #read in metadata to find latest time slot. Then assign this to FIRSTDAY
        METADATA_PATH = '{}metadata.pickle'.format(OUTFILESPATH.\
                  replace(os.path.basename(OUTFILESPATH), ""))
      
        metadata = pscrosscorr.load_pickled_xcorr(METADATA_PATH)
        #print "metadata: ", metadata[-5:]
        #re-assign FIRSTDAY variable to where the data was cut off
        #del metadata[-1]
        FIRSTDAY = metadata[len(metadata) - 1] #+  \
        #dt.timedelta(minutes=XCORR_INTERVAL)
        
    # FIND RESTART DATE FROM PARTIAL PICKLE FILE, NOT THE METADATA PICKLE
    
    
    # ============
    # Main program
    # ============
    
    # Reading inventories in dataless seed and/or StationXML files
    dataless_inventories = []
    if USE_DATALESSPAZ:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dataless_inventories = psstationSQL.get_dataless_inventories(DATALESS_DIR,
                                                                      verbose=False)
    
    xml_inventories = []
    if USE_STATIONXML:
        xml_inventories = psstationSQL.get_stationxml_inventories(STATIONXML_DIR,
                                                               verbose=False)
    
    # Getting list of stations
    #stations, subdir_len = psstation.get_stations(mseed_dir=MSEED_DIR,
    #                                  xml_inventories=xml_inventories,
    #                                  dataless_inventories=dataless_inventories,
    #                                  startday=FIRSTDAY,
    #                                  endday=LASTDAY,
    #                                  verbose=False)
    
    #connect SQL database
    SQL_db = os.path.join(DATABASE_DIR, 'timeline.db')
    
    stations, subdir_len = psstationSQL.get_stationsSQL(SQL_db, 
                           xml_inventories=xml_inventories,
                           dataless_inventories=dataless_inventories,
                           startday=FIRSTDAY,
                           endday=LASTDAY,
                           verbose=False)

    print "stations: ", stations

    stat_coords = np.asarray([station.coord for station in stations])

    
    DECLUSTER = False
    
    #if DECLUSTER: 
    #    stat_coords = np.asarray([station.coord for station in stations])
    #    COORDS = Coordinates(input_list=stat_coords)
    #    declustered_coords = COORDS.decluster(degree_dist=0.1)
    
    #    stations = [station for station in stations if 
    #                station.coord in declustered_coords]      


    # Loop on time interval
     #number of time steps
    N = int(((LASTDAY - FIRSTDAY).days + 1)*60*24 / XCORR_INTERVAL)
    
    dates = [FIRSTDAY + dt.timedelta(minutes=i) for i in \
             [j*XCORR_INTERVAL for j in range(N)]]
    
    if RANDOM_STACK:
        #print dates        
        # randomly shuffle dates 
        #dates = np.array(dates)
        np.random.shuffle(dates)
        #print dates
        
    
    #begin = raw_input("\nPress enter to begin the program ")
    
    # initialise preprocess class: METHOD - Bensen et al. (2007)
    Preprocess = pspreprocess.Preprocess(FREQMIN, FREQMAX,
                                         FREQMIN_EARTHQUAKE,
                                         FREQMAX_EARTHQUAKE,
                                         CORNERS, ZEROPHASE,
                                         PERIOD_RESAMPLE,
                                         WINDOW_TIME,
                                         WINDOW_FREQ,
                                         ONEBIT_NORM)
                          
    #loop on time-series. Date now represents XCORR_INTERVAL long time intervals
    counter = 0 
    
    
    for date in dates:
        print date
        loop_time0 = dt.datetime.now()

        print "\nProcessing data for date {} with a {} minute cross-correlation\
 time-interval between times: {} and {}".format(date.date(), \
        int(XCORR_INTERVAL)  , date.time(), \
        (date + dt.timedelta(minutes=XCORR_INTERVAL)).time())
        
        iterate_stations = sorted(sta for sta in stations)
        
        # =====================================================================
        # check iterate stations have a file in the SQL database
        # =====================================================================
        
        # connect the database
        conn = lite.connect(SQL_db)
        # create cursor object
        c = conn.cursor()
        
        # convert to UTC timestamp to search in SQL database
        search_start = (date - dt.timedelta(minutes=1) - epoch).total_seconds()
        search_end =  (date + dt.timedelta(minutes=XCORR_INTERVAL+1) - epoch).total_seconds()     
        
        # check if files have data within the time frame search_end-search_start
        populated_stations = c.execute('SELECT station FROM file_extrema WHERE \
starttime <= ? AND endtime >= ?', (search_start, search_end))

        populated_stations = list(it.chain(*list(populated_stations.fetchall())))
        
        # filter stations with no data for the given time period of this loop! 
        for stat in iterate_stations:
            stat_code = unicode('{}.{}.{}'.format(stat.network, 
                                                  stat.name, 
                                                  stat.channel))
            
            if stat_code not in populated_stations:
                iterate_stations.remove(stat)


        # close timeline.db database
        conn.close()

        iterate_stations = iterate_stations[1:]
        # =====================================================================
        # =====================================================================            

        # subset if stations (if provided)
        if CROSSCORR_STATIONS_SUBSET:
            iterate_stations = [sta for sta in iterate_stations
                              if sta.name in CROSSCORR_STATIONS_SUBSET]
                                  
        # =================
        # processing traces
        # =================                          
    
        # =============================================================
        # preparing functions that get one merged trace per station
        # and pre-process trace, ready to be parallelized (if required)
        # =============================================================

        def get_merged_trace(station):
            """
            Preparing func that returns one trace from selected station,
            at current date. Function is ready to be parallelized.
            """

            try:
                trace = Preprocess.get_merged_trace(station=station,
                                                     date=date,
                                                     xcorr_interval=XCORR_INTERVAL,
                                                     skiplocs=CROSSCORR_SKIPLOCS,
                                                     minfill=MINFILL)
                                                     
                
                #plt.figure()
                #plt.plot(trace.data)
                #plt.show()
                #plt.clf()
                
                
                if total_verbose:
                    msg = 'merged'
                    print '{}.{} [{}] '.format(trace.stats.network, 
                                            trace.stats.station, 
                                            msg),
                errmsg = None
                    
            except pserrors.CannotPreprocess as err:
                # cannot preprocess if no trace or daily fill < *minfill*
                                
                trace = None
                errmsg = '{}: skipping'.format(err)
            except Exception as err:
                # unhandled exception!
                trace = None
                errmsg = 'Unhandled error: {}'.format(err)
    
            if errmsg:
                # printing error message
                if total_verbose:
                    print '{}.{} [{}] '.format(station.network, 
                                               station.name, errmsg),
    
            return trace
    
        def preprocessed_trace((trace, response)):
            """
            Preparing func that returns processed trace: processing includes
            removal of instrumental response, band-pass filtering, demeaning,
            detrending, downsampling, time-normalization and spectral whitening
            (see pscrosscorr.preprocess_trace()'s doc)
    
            Function is ready to be parallelized.
            """
            


            #if not trace or response is False:
            #    return
            
    	    if not type(trace) is Trace:
		print trace
		#plt.figure()
		#plt.plot(trace)
		#plt.show()
		#quit()
		return


        #try:

            Preprocess.preprocess_trace(trace=trace, paz=response, 
                                        verbose=False)
            msg = 'ok'
            if total_verbose:
                print '{}.{} [{}] '.format(trace.stats.network, 
                                            trace.stats.station, 
                                            msg),
#    
#            except pserrors.CannotPreprocess as err:
#                # cannot preprocess if no instrument response was found,
#                # trace data are not consistent etc. (see function's doc)
#                trace = None
#                print(err)
#                print 'skipping'
#    
#            except Exception as err:
#                # unhandled exception!
#                trace = None
#                print(err)
                print 'skipping'
    
            # printing output (error or ok) message
    
            # although processing is performed in-place, trace is returned
            # in order to get it back after multi-processing
            return trace
    
        # ====================================
        # getting one merged trace per station
        # ====================================
        merge_t0 = dt.datetime.now()
        print '\nMerging traces ... '
        if MULTIPROCESSING['merge trace']:
            # multiprocessing turned on: one process per station
            pool = mp.Pool(None)
            traces = pool.map(get_merged_trace, iterate_stations)
            pool.close()
            pool.join()
        else:
            # multiprocessing turned off: processing stations one after another
            traces = [get_merged_trace(s) for s in iterate_stations]
    
        # =====================================================
        # getting or attaching instrumental response
        # (parallelization is difficult because of inventories)
        # =====================================================
        
        #print "traces1: ", traces
        responses = []
        for tr in traces:
            
            if not tr:
                responses.append(None)
                continue
    
            # responses elements can be (1) dict of PAZ if response found in
            # dataless inventory, (2) None if response found in StationXML
            # inventory (directly attached to trace) or (3) False if no
            # response found
            if RESP_REMOVE:  
                try:
                    response = Preprocess.get_or_attach_response(
                        trace=tr,
                        dataless_inventories=dataless_inventories,
                        xml_inventories=xml_inventories)
                    errmsg = None
                except pserrors.CannotPreprocess as err:
                    # response not found
                    response = False
                    errmsg = '{}: skipping'.format(err)
                except Exception as err:
                    # unhandled exception!
                    response = False
                    errmsg = 'Unhandled error: {}'.format(err)
    
                responses.append(response)
                if errmsg:
                    # printing error message
                    if total_verbose:
                        print '{}.{} [{}] '.format(tr.stats.network, 
                                                   tr.stats.station, 
                                                   errmsg),
            else:
                responses.append(None)
                                           
        #print "traces2: ", traces
        print '\nTraces merged and responses removed in {:.1f} seconds'\
        .format((dt.datetime.now() - merge_t0).total_seconds()) 

        # =================
        # processing traces
        # =================
        print '\nPre-processing traces ... '

        t0 = dt.datetime.now()
        # must have more than one trace for cross-correlations!
        #traces = np.array(traces)
        #traces_check = traces[traces != np.array(None)]
        #print "traces3: ", traces


        if MULTIPROCESSING['process trace']:
                # multiprocessing turned on: one process per station
                pool = mp.Pool(NB_PROCESSES)
                traces = pool.map(preprocessed_trace, zip(traces, responses))
                pool.close()
                pool.join()
            
        else:
                # multiprocessing turned off: processing stations one after another
                try:
                    traces = map(preprocessed_trace, zip(traces, responses))
                except:
                    continue
        # setting up dict of current date's traces, {station: trace}
        tracedict = {s.name: trace for s, trace in zip(iterate_stations, 
                                                       traces) if trace}
    
        delta = (dt.datetime.now() - t0).total_seconds()
        print "\nProcessed traces in {:.1f} seconds".format(delta)

        # create tmp folder for tracedict
        #if not os.path.exists('tmp'): os.makedirs('tmp')   
        #dump the time interval's pre-processed items in tracedict to a pickle
        #with open('tmp/preprocessed_tracedict.pickle', 'wb') as f:
        #    print "\nExporting pre-processed traces of time-series to: " + f.name
        #    pickle.dump(tracedict, f, protocol=2)
        
                                       
        # import tracedict from output pickle produced with preprocess_total
        #tracedict_pickle = 'tmp/preprocessed_tracedict.pickle'
        #f = open(name=tracedict_pickle, mode='rb')
        #tracedict = pickle.load(f)
        #f.close()
        
        # remove preprocessed tracedict pickle file
        #if os.path.isfile(tracedict_pickle): os.remove(tracedict_pickle)
            
            
        # ======================================================
        # stacking cross-correlations of the current time-series
        # ======================================================
    
        
        #if len(tracedict) < 2:
        #    print "No cross-correlation for this interval"
        #    continue
    
        t0 = dt.datetime.now()
        xcorrdict = {}
        if MULTIPROCESSING['cross-corr']:
            # if multiprocessing is turned on, we pre-calculate cross-correlation
            # arrays between pairs of stations (one process per pair) and feed
            # them to xc.add() (which won't have to recalculate them)
            print "\nProcessing cross-correlations ..."
    
            def xcorr_func(pair):
                """
                Preparing func that returns cross-correlation array
                beween two traces
                """
                (s1, tr1), (s2, tr2) = pair
               
                print '{}-{} '.format(s1, s2),
                shift = int(CROSSCORR_TMAX / PERIOD_RESAMPLE)
                xcorr = obspy.signal.cross_correlation.xcorr(
                    tr1, tr2, shift_len=shift, full_xcorr=True)[2]
                    
                #plt.figure()
                #plt.title("xcorr 1 ")
                #plt.plot(xcorr)
                #plt.show()
                #plt.clf()
                return xcorr
    

            pairs = list(it.combinations(sorted(tracedict.items()), 2))
            pool = mp.Pool(NB_PROCESSES)
            xcorrs = pool.map(xcorr_func, pairs)
            pool.close()
            pool.join()
            xcorrdict = {(s1, s2): xcorr for ((s1, _), (s2, _)), 
                         xcorr in zip(pairs, xcorrs)}
            print 

        print "Stacking cross-correlations"
        xc.add(tracedict=tracedict,
               stations=stations,
               xcorr_tmax=CROSSCORR_TMAX,
               xcorrdict=xcorrdict,
               date=date,
               verbose=not MULTIPROCESSING['cross-corr'])
        
        pairs = list(it.combinations(sorted(tracedict.items()), 2))

        # calculate max snr for snr weighted stack! 
#        for pair in pairs: 
#            (s1, tr1), (s2, tr2) = pair
#            s1, s2 = str(s1), str(s2)
#            snr_list = xc[s1][s2].SNR_list
#            max_snr = np.max(snr_list)

#            snr_stack = xc[s1][s2].SNR_stack
            
#            snr_wstack = np.zeros_like(snr_stack[0])
             
#            for xcorr, snr in zip(snr_stack, snr_list): 
#                snr_wstack += xcorr * snr / max_snr
                
            # assign final snr weighted stack xcorr green's function to SNR_stack
#            xc[s1][s2].SNR_stack = snr_wstack
#            if s1 in xc.keys():
#                if s2 in xc[s1].keys():
#                    pws = xc[s1][s2].pws
#                    plt.figure()
#                    plt.plot(pws)
#                    plt.show()
#                    plt.clf()

            
    
    #==============================================================================    
        delta = (dt.datetime.now() - t0).total_seconds()


        #print "\nCalculated and stacked {} cross-correlations in \
#{:.1f} seconds".format(len(xcorrs), delta)
    
        loop_delta = (dt.datetime.now() - loop_time0).total_seconds()

        print "\nIn total, the previous loop had a processing time of: \
{:.1f} seconds".format(loop_delta)
             
            #there may be an overlap in metadata times. Maybe add xcorr interval to this?
        #have the program restart from the beginning of each stack! not each month
        #this gives more redundancy!
        #(allows to restart after a crash from that date)
        
        
        # save partial pickle file only if timeseries loop is large enough
        #if len(metadata) >= 2:
        #    print "Time since last save: ", abs(date - metadata[-1] + dt.timedelta(minutes=XCORR_INTERVAL)) 
            #if len(metadata) % 10 == 0:# (date - metadata[-1]) >= dt.timedelta(hours=1):
        #    with open(u'{}.part.pickle'.format(OUTFILESPATH), 'wb') as f:
        #        print "\nExporting cross-correlations calculated until now." 
        #        pickle.dump(xc, f, protocol=2)
        #    metadata.append(date)
                
        #elif len(metadata) == 0:
        #    metadata.append(date)

    
            
        #also create a metadata dump file for use only if the program needs to be restarted
        #use replace() to get rid of basename to create file named metadata.pickle in 
        #correct path
            
        #with open(u'{}.part.pickle'.format(OUTFILESPATH), 'wb') as f:
        #    print "\nExporting cross-correlations calculated until now." 
        #    pickle.dump(xc, f, protocol=2)    
        with open(METADATA_PATH, 'wb') as f:
            print "\nExporting re-start metadata of time-series calculated until \
now."
            pickle.dump(metadata, f, protocol=2)
    
    # export xcorr pairs
    if not xc.pairs():
        print "No cross-correlation could be calculated: nothing to export!"
    else:
    	    # exporting to binary and ascii files
        xc.export(outprefix=OUTFILESPATH, stations=stations, verbose=True)
    

    if not xc.pairs():
        print "No cross-correlation could be calculated: nothing to export!"
    else:
    	    # exporting to binary and ascii files
        xc.export(outprefix=OUTFILESPATH+'.part', stations=stations, verbose=True)
        # exporting to png file
    print "Exporting cross-correlations to file: {}.png".format(OUTFILESPATH)
        # optimizing time-scale: max time = max distance / vmin (vmin = 2.5 km/s)
    maxdist = max([xc[s1][s2].dist() for s1, s2 in xc.pairs()])
    maxt = min(CROSSCORR_TMAX, maxdist / 2.5)
        

    if PLOT_DISTANCE:
            #plot distance plot of cross-correlations
        xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
                    outfile=os.path.join(OUTFOLDERS, OUTFILESNAME)\
                    + '.png', showplot=False, stack_type='linear')
                    
        #xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
        #            outfile=os.path.join(OUTFOLDERS, OUTFILESNAME)\
        #            + '.png', showplot=False, stack_type='PWS')
                    
        #xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
        #            outfile=os.path.join(OUTFOLDERS, OUTFILESNAME)\
        #            + '.png', showplot=False, stack_type='SNR')
                    
                    
        #    xc.plot(plot_type='distance', xlim=(-maxt, maxt), 
        #            outfile=os.path.join(OUTFOLDERS, OUTFILESNAME)\
        #            + '.png', showplot=False, stack_type='combined')                    
                    
                    
    if PLOT_CLASSIC:
            #plot individual cross-correlations
        xc.plot(plot_type='classic', xlim=(-maxt, maxt), 
                    outfile=OUTFOLDERS, showplot=False)

        
        #xc.plot_SNR(plot_type='individual', outfile=OUT_SNR)

    # removing file containing periodical exports of cross-corrs
    # only do this if the full file exists!     
    try:
        os.remove(u'{}.part.pickle'.format(OUTFILESPATH))
    except:
        pass


    #remove_config(config_file)
    

total_delta = (dt.datetime.now() - total_time0).total_seconds()


print "Calculated every xcorr in time-series in in \
{:.1f} seconds".format(total_delta)
