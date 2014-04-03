import numpy as np

def W3FB05(XI,XJ,XMESHL=60,ORIENT=40):
    """
    C$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
    C
    C SUBPROGRAM: W3FB05         GRID COORDINATES TO LATITUDE, LONGITUDE
    C   AUTHOR: JONES,R.E.       ORG: W345       DATE: 86-07-17
    C
    C ABSTRACT: CONVERTS THE COORDINATES OF A LOCATION FROM THE GRID(I,J)
    C   COORDINATE SYSTEM OVERLAID ON THE POLAR STEREOGRAPHIC MAP PROJEC-
    C   TION TRUE AT 60 DEGREES N OR S LATITUDE TO THE NATURAL COORDINATE
    C   SYSTEM OF LATITUDE/LONGITUDE ON THE EARTH. W3FB05 IS THE REVERSE
    C   OF W3FB04.
    C
    C PROGRAM HISTORY LOG:
    C   86-07-17  R.E.JONES
    C   89-11-01  R.E.JONES   CHANGE TO CRAY CFT77 FORTRAN
    C
    C USAGE:  CALL W3FB05 (XI, XJ, XMESHL, ORIENT, ALAT, ALONG)
    C
    C   INPUT VARIABLES:
    C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
    C     ------ --------- -----------------------------------------------
    C     XI     ARG LIST  I OF THE POINT RELATIVE TO THE NORTH OR S. POLE
    C     XJ     ARG LIST  J OF THE POINT RELATIVE TO THE NORTH OR S. POLE
    C     XMESHL ARG LIST  MESH LENGTH OF GRID IN KM AT 60 DEGREES(<0 IF SH)
    C                   (190.5 LFM GRID, 381.0 NH PE GRID,-381.0 SH PE GRID)
    C     ORIENT ARG LIST  ORIENTATION WEST LONGITUDE OF THE GRID
    C                    (105.0 LFM GRID, 80.0 NH PE GRID, 260.0 SH PE GRID)
    C
    C   OUTPUT VARIABLES:
    C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
    C     ------ --------- -----------------------------------------------
    C     ALAT   ARG LIST  LATITUDE IN DEGREES  (<0 IF SH)
    C     ALONG  ARG LIST  WEST LONGITUDE IN DEGREES
    C
    C   SUBPROGRAMS CALLED:
    C     NAMES                                                   LIBRARY
    C     ------------------------------------------------------- --------
    C     ASIN   ATAN2                                            SYSLIB
    C
    C   REMARKS: ALL PARAMETERS IN THE CALLING STATEMENT MUST BE
    C     REAL. THE RANGE OF ALLOWABLE LATITUDES IS FROM A POLE TO
    C     30 DEGREES INTO THE OPPOSITE HEMISPHERE.
    C     THE GRID USED IN THIS SUBROUTINE HAS ITS ORIGIN (I=0,J=0)
    C     AT THE POLE, SO IF THE USER'S GRID HAS ITS ORIGIN AT A POINT
    C     OTHER THAN A POLE, A TRANSLATION IS REQUIRED TO GET I AND J FOR
    C     INPUT INTO W3FB05. THE SUBROUTINE GRID IS ORIENTED SO THAT
    C     GRIDLINES OF I=CONSTANT ARE PARALLEL TO A WEST LONGITUDE SUP-
    C     PLIED BY THE USER. THE EARTH'S RADIUS IS TAKEN TO BE 6371.2 KM.
    C
    C   WARNING: THIS CODE WILL NOT VECTORIZE, IT IS NORMALY USED IN A
    C            DOUBLE DO LOOP WITH W3FT01, W3FT00, ETC. TO VECTORIZE IT,
    C            PUT IT IN LINE, PUT W3FT01, W3FT00, ETC. IN LINE.
    C
    C   LANGUAGE: CRAY CFT77 FORTRAN
    C   MACHINE:  CRAY Y-MP8/832
    C
    C$$$
    converted from Fortran to Python by Timothy W. Hilton, 1 April 2014

    fortran source from
    http://www.nco.ncep.noaa.gov/pmb/docs/libs/w3lib/w3fb05.html,
    http://www.nco.ncep.noaa.gov/pmb/docs/libs/w3lib/w3lib.tar,
    accessed 1 April 2014
    """
    DEGPRD = 57.2957795
    EARTHR = 6371.2

    GI2   = ((1.86603 * EARTHR) / (XMESHL))**2
    R2    = XI * XI + XJ * XJ

    if (np.abs(R2) < 1e-6):
        ALONG = 0.0
        ALAT  = 90.0
        if (XMESHL < 0.0):
            ALAT = -ALAT
        return(ALONG, ALAT)
    else:
        ALAT  = np.arcsin((GI2 - R2) / (GI2 + R2)) * DEGPRD
        ANGLE = DEGPRD * np.arctan2(XJ,XI)
        if (ANGLE < 0.0):
             ANGLE = ANGLE + 360.0


    if (XMESHL > 0.0):
        ALONG = 270.0 + ORIENT - ANGLE
    else:
        ALONG = ANGLE + ORIENT - 270.0
        ALAT  = -(ALAT)

    if (ALONG < 0.0):
        ALONG = ALONG + 360.0 
    if (ALONG >= 360.0):
        ALONG = ALONG - 360.0

    return(ALONG, ALAT)


if __name__ == "__main__":
    import mpl_toolkits.basemap as basemap
    import matplotlib.pyplot as plt

    # make a 124 x 124 grid of indices, 1 to 124
    xi, xj = np.meshgrid(np.arange(124) + 1, np.arange(124) + 1)
    # calculate lon, lat for each grid cell
    gridsize = 60
    reflon = 40
    lon = np.zeros(xi.shape)
    lat = np.zeros(xi.shape)
    for i in np.arange(xi.shape[0]):
        for j in np.arange(xi.shape[1]):
            lon[i,j], lat[i,j] = W3FB05(xi[i,j], xj[i,j], gridsize, reflon)

    #=====
    #draw the lon, lat pairs on a map
    #=====

    #set up the map, draw labels, etc.
    mapwidth = 12.0e6  # not sure of units for width/height
    mapheight = 9.5e6
    R_EARTH = 6371007.181000  #Earth radius in meters
    m = basemap.Basemap(width=mapwidth,
                        height=mapheight,
                        projection='aeqd',
                        lat_0=54,
                        lon_0=-105,
                        resolution='l',
                        area_thresh=1000, #show features larger than 1000 km
                        rsphere=R_EARTH,
                        fix_aspect=True)
    m.drawcoastlines()
    m.drawmeridians(meridians=range(0, -180, -20),
                    labels=(0, 0, 0, 1))  # labels on bottom
    m.drawparallels(circles=range(0, 90, 20),
                    labels=(1, 0, 0, 0))  # labels on left
    #draw the STEM grid
    m.plot(lon,
           lat,
           latlon=True,
           marker='.',
           markersize=0.2,
           color='blue')
    plt.title('lon, lat from W3FB05')

    plt.show()

