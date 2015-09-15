import pandas as pd
import numpy as np
import MySQLdb
import socket


def select_kics(catfile='villanova-db.csv', pmin=0.0, pmax=None):
    """
    Return KIC IDs based on system parameters.

    Inputs
    ----------
    catfile : string
        Name of catalog file. (Default: villanova-db.csv)
    pmin : float, optional
        Minimum orbital period in days. (Default: 0.0)
    pmax : float, optional
        maximum orbital period in days.

    Returns
    -------
    kics : array-like
        KIC IDs matching search criteria.

    """
    # Construct the query string.
    qstring = 'period > %s & ' % str(pmin)
    if pmax:
        qstring += 'period < %s & ' % str(pmax)

    # Load the Villanova catalog, find the matching KIC IDs.
    df = pd.read_csv(catfile)
    kics = df.query(qstring[:-2])['KIC'].values.astype(int).astype(str)

    print 'Found %d systems in catalog meeting criteria.' % kics.size
    return kics.astype(int)


def __dbconnect(host='tddb.astro.washington.edu', user='tddb', password='tddb',
                db='Kepler'):
    """
    Log into a database using MySQLdb. Written by Ethan Kruse.

    Parameters
    ----------
    host : string, optional
        Default tddb.astro.washington.edu
    user : string, optional
        Default eakruse
    password : string, optional
        Default tddb
    db : string, optional
        Default Kepler

    Returns
    -------
    dbconnect : Connect
        MySQLdb connector.

    """
    return MySQLdb.connect(host=host, user=user, passwd=password, db=db,
                           connect_timeout=0)


def loadlc_db(kic, usepdc=True, lc=True, **kwargs):
    """
    Load Kepler data from the local tddb database. Written by Ethan Kruse.

    Can pass optional MySQLdb keyword arguments for logging into the database
    (host, user, passwd, db). Those default values should work though.

    Parameters
    ----------
    kic : int
        Kepler Input Catalog number for the target.
    usepdc : bool, optional
        Defaults to True. If True, use the PDCSAP data
        instead of the raw SAP.
    lc : bool, optional
        Whether to select long or short cadence. Defaults to True, or LC data.

    Returns
    -------
    time : ndarray
        Kepler times of center of exposure
    flux : ndarray
        Kepler fluxes normalized for each quarter
    fluxerr : ndarray
        Kepler flux errors for each exposure
    cadence : ndarray
        Cadence number
    quarter : ndarray
        Kepler quarter
    quality : ndarray
        Kepler data quality flag

    """
    tablename = 'source'
    if lc:
        lcflag = "LCFLAG > 0"
    else:
        lcflag = "LCFLAG = 0"

    if usepdc:
        fluxstr = "pdcsap_flux, pdcsap_flux_err "
    else:
        fluxstr = "sap_flux, sap_flux_err "

    hostname = socket.gethostname()
    if 'astro.washington.edu' in hostname:
        ct = 0
        gotit = False
        # try multiple times in case of sporadic database timeouts
        while ct < 5 and not gotit:
            try:
                db = __dbconnect(**kwargs)
                cursor = db.cursor()

                toex = "SELECT cadenceno, quarter, sap_quality, time, {0} FROM" \
                       " {2} WHERE keplerid = %s AND {1};"\
                    .format(fluxstr, lcflag, tablename)

                cursor.execute(toex, (int(kic),))
                results = cursor.fetchall()
                cadence = np.array([x[0] for x in results], dtype=np.int32)
                quarter = np.array([x[1] for x in results], dtype=np.int32)
                quality = np.array([x[2] for x in results], dtype=np.int32)
                time = np.array([x[3] for x in results], dtype=np.float64)
                flux = np.array([x[4] for x in results], dtype=np.float32)
                fluxerr = np.array([x[5] for x in results], dtype=np.float32)
                cursor.close()
                db.close()
                # for some reason some results are coming back with arrays
                # of length 0.
                if len(time) > 0:
                    gotit = True
                ct += 1
            except MySQLdb.OperationalError:
                print "mysqldb connection failed on attempt {0} of {1}.\n" \
                      "Trying again.".format(ct + 1, 5)
                ct += 1
    else:
        import paramiko

        toex = "SELECT cadenceno, quarter, sap_quality, time, {0} FROM {3} " \
               "WHERE keplerid = {2} AND {1};"\
            .format(fluxstr, lcflag, int(kic), tablename)
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        # no password because my SSH keygen doesn't have a password
        ssh.connect('hail.astro.washington.edu', username='eakruse')
        stdin, stdout, stderr = ssh.exec_command(
            'mysql -h tddb.astro.washington.edu -D Kepler -u eakruse '
            '--password=tddbKepler -e "{0}"'.format(toex))
        results = stdout.read().splitlines()
        results = results[1:]

        cadence = np.array([int(x.split('\t')[0]) for x in results],
                           dtype=np.int32)
        quarter = np.array([int(x.split('\t')[1]) for x in results],
                           dtype=np.int32)
        quality = np.array([int(x.split('\t')[2]) for x in results],
                           dtype=np.int32)
        time = np.array([float(x.split('\t')[3]) for x in results],
                        dtype=np.float64)
        flux = np.array([float(x.split('\t')[4]) for x in results],
                        dtype=np.float32)
        fluxerr = np.array([float(x.split('\t')[5]) for x in results],
                           dtype=np.float32)
        ssh.close()

    # guarantee the light curve is in sequential order
    # %timeit says that doing the ordering in python is faster than including
    # an 'ORDER BY time' flag in the mysql search. I have no idea why, but
    # I'll keep doing the ordering here.
    order = np.argsort(time)
    time = time[order]
    flux = flux[order]
    fluxerr = fluxerr[order]
    quality = quality[order]
    cadence = cadence[order]
    quarter = quarter[order]

    # go from raw CCD counts to normalized fluxes per quarter
    uquarts = np.unique(quarter)
    for ii in uquarts:
        val = np.where(quarter == ii)[0]
        fluxerr[val] /= np.median(flux[val])
        flux[val] /= np.median(flux[val])

    if flux.size == 0:
        print 'No light curves found!'

    return time, flux, fluxerr, cadence, quarter, quality
