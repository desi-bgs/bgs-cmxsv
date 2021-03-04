'''

module for quickly calculating observing conditions for kitt peak 

'''
import ephem
import astropy.units as u
import astropy.time as atime

mayall = ephem.Observer()
mayall.lat = (31.963972222 * u.deg).to(u.rad).value
mayall.lon = (-111.599336111 * u.deg).to(u.rad).value
mayall.elevation = 2120. #m 
mayall.pressure = 1e3 * (78318 * u.Pascal).to(u.bar).value # mbar 
mayall.temp = 5. #C 


def get_sun(mjd): 
    ''' return sun ra, dec, and altitude given mjd for Mayall

    Parameters
    ----------
    mjd : float
        mjd date

    Returns
    -------
    ra : float
        RA of the sun at mjd in degrees

    dec : float
        Dec of sun at mjd in degrees

    alt : float
        altitude of the sun  in degrees
    '''
    obs = mayall.copy() 
    obs.date = atime.Time(mjd, format='mjd').iso

    sun = ephem.Sun()
    sun.compute(obs)
    
    ra = (sun.ra * u.radian).to(u.degree).value 
    dec = (sun.dec * u.radian).to(u.degree).value
    alt = (sun.alt * u.radian).to(u.degree).value
    return ra, dec, alt


def get_moon(mjd): 
    ''' return sun ra, dec, and altitude given mjd for Mayall

    Parameters
    ----------
    mjd : float
        mjd date

    Returns
    -------
    ra : float
        RA of the sun at mjd in degrees

    dec : float
        Dec of sun at mjd in degrees

    moon_phase : float 
        moon illumination fraction

    alt : float
        altitude of the sun  in degrees
    '''
    obs = mayall.copy() 
    obs.date = atime.Time(mjd, format='mjd').iso

    moon = ephem.Moon()
    moon.compute(obs)
    
    ra = (moon.ra * u.radian).to(u.degree).value 
    dec = (moon.dec * u.radian).to(u.degree).value
    alt = (moon.alt * u.radian).to(u.degree).value
    return ra, dec, moon.moon_phase, alt


def ang_separation(ra1, dec1, ra2, dec2): 
    ''' angular separation between (ra1, dec1) and (ra2, dec2) 

    Parameters
    ----------
    ra1 : float
        RA in degrees

    dec1 : float 
        Dec in degrees

    ra2 : float
        RA in degrees

    dec2 : float 
        Dec in degrees

    Returns
    -------
    sep : float
        separation in degrees
    '''
    ra1, dec1 = np.deg2rad(ra1), np.deg2rad(dec1)
    ra2, dec2 = np.deg2rad(ra2), np.deg2rad(dec2)

    havRA12 = 0.5 * (1 - np.cos(ra2 - ra1))
    havDEC12 = 0.5 * (1 - np.cos(dec2 - dec1))
    havPHI = havDEC12 + np.cos(dec1) * np.cos(dec2) * havRA12
    return np.rad2deg(np.arccos(np.clip(1 - 2 * havPHI, -1, +1)))
