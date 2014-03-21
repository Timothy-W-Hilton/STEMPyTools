from scipy.io import netcdf
import re
import math, os
import numpy as np
from datetime import datetime
import pandas as pd
import csv
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
from scipy.spatial import cKDTree
from brewer2mpl import qualitative as brewer_qualitative
import STEM_parsers

class NOAA_OCS(object):
    """
    Container class to work with NOAA airborne OCS observations
    """
    def __init__(self,
                 obs=None):
        self.obs = obs
        self.obs_color = brewer_qualitative.Dark2['max'].mpl_colors[0]
        self.grid_color = brewer_qualitative.Dark2['max'].mpl_colors[1]

    @classmethod
    def parse_file(cls, fname):
        """
        parses a NOAA airborne OCS measurement text file.
        PARAMETERS
        ----------
        fname: full path the NOAA OCS file to be parsed
        """

        # the NOAA files contain varying number of header lines.  It
        # looks like the column names are always the last commented
        # header line.  Now get the column names and the number of 
       # header lines.
        f = open(fname, 'r')
        ln = "#"  #initialize to empty comment
        hdr_count = 0
        while(ln[0] is '#'):
            last_ln = ln
            ln = f.readline()
            hdr_count += 1
        names = re.sub('.*data_fields:\s+', '', last_ln).rstrip().split()
        # parse the data to a DataFrame
        df = pd.read_csv(fname,
                         comment='#',
                         delim_whitespace=True,
                         skiprows=hdr_count-1,
                         names=names)
        # convert the times to datetime objects
        df['date'] = df.apply(
            lambda row: datetime(row['sample_year'],
                                 row['sample_month'],
                                 row['sample_day'],
                                 row['sample_hour'],
                                 row['sample_minute'],
                                 row['sample_seconds']),
            axis=1)

        return(cls(obs=df))

    def get_stem_xy(self, stem_lon, stem_lat):
        """
        Get STEM grid (x, y) values for the object's observations.
        PARAMETERS
        ----------
        stem_lon: longitudes of STEM grid cell centers
        stem_lat: latitudes of STEM grid cell centers
        """
        self.obs['x_stem'], self.obs['y_stem'] = get_stem_xy(
            self.obs.sample_longitude,
            self.obs.sample_latitude,
            stem_lon,
            stem_lat)
        self.obs['lon_stem'] = stem_lon[self.obs.x_stem, self.obs.y_stem]
        self.obs['lat_stem'] = stem_lat[self.obs.x_stem, self.obs.y_stem]

    def calcT(start, obsT):
        """
        Calculates the stem t. Start is the start of the stem simulation
        and obsT is the observation time. Both are given as [hour, day,
        month, year].
        """
        SECONDS_PER_HOUR = 60*60
        # convert start, obsT to datetime.datetime
        start = datetime(start[3], start[2], start[1], start[0], 0, 0)
        obsTt = datetime(obsT[3], obsT[2], obsT[1], obsT[0], 0, 0)

        dt = (obsT - start).total_seconds() / SECONDS_PER_HOUR
        return(dt)

    def plot_stem_obs_matches(self):
        """
        plot observations (^) and nearest STEM grid cells (o) on a
        map, with a line connecting each observation to its assigned
        grid cell
        """
        lats = np.array((self.obs.sample_latitude.values,
                         self.obs.lat_stem.values))
        lons = np.array((self.obs.sample_longitude.values,
                         self.obs.lon_stem.values))
        lat_min = lats.min()
        lat_max = lats.max()
        lon_min = lons.min()
        lon_max = lons.max()
        lat_pad = (lat_max - lat_min) * 0.1
        lon_pad = (lon_max - lon_min) * 0.1

        # lon, lat for SW corner of map
        SW_crnr = (lon_min - lon_pad, lat_min - lat_pad) 
        # lon, lat for NE corner of map
        NE_crnr = (lon_max + lon_pad, lat_max + lat_pad) 
        m = basemap.Basemap(llcrnrlon=SW_crnr[0], llcrnrlat=SW_crnr[1],
                            urcrnrlon=NE_crnr[0], urcrnrlat=NE_crnr[1])

        m.drawmeridians(np.linspace(lon_max, lon_min, num=3),
        			labels=[1,1,1,1])
        m.drawparallels(np.linspace(lat_max, lat_min, num=3),
        			labels=[1,1,1,1])

        obs_mrks = m.scatter(self.obs.sample_longitude,
                             self.obs.sample_latitude,
                             c=self.obs_color,
                             label='observation',
                             zorder=2)
        grid_mrks = m.scatter(self.obs['lon_stem'].values,
                              self.obs['lat_stem'].values,
                              color=self.grid_color,
                              marker='^',
                              label='STEM grid cell',
                              zorder=3)
        n_obs = self.obs.shape[0]
        for i in range(n_obs):
            ln = m.plot((self.obs['lon_stem'].values[i],
                         self.obs['sample_longitude'].values[i]),
                        (self.obs['lat_stem'].values[i],
                         self.obs['sample_latitude'].values[i]),
                        color='black',
                        zorder=1)
        plt.legend(numpoints=1,
                   scatterpoints=1,
                   loc='best')
        return(m)
        
    def plot_locations(self):
        m = init_NA_map()
        mrks = m.scatter(
            self.obs.sample_longitude,
            self.obs.sample_latitude,
            c=self.obs_color,
            label='observation')
        return(m, mrks)

def init_NA_map():
    """
    plot data locations on a map of North America
    """
    plt.figure(figsize=(8,8))
    #initialize a map of North America
    SW_crnr = (-170, 10) # lon, lat for SW corner of map
    NE_crnr = (-50, 80) # lon, lat for NE corner of map
    m = basemap.Basemap(llcrnrlon=SW_crnr[0], llcrnrlat=SW_crnr[1],
                        urcrnrlon=NE_crnr[0], urcrnrlat=NE_crnr[1])
    m.drawmeridians(np.linspace(SW_crnr[0], NE_crnr[0], num=3),
			labels=[1,1,1,1])
    m.drawparallels(np.linspace(SW_crnr[1], NE_crnr[1], num=3),
			labels=[1,1,1,1])
    m.drawcoastlines()
    return(m)
    
def get_height_levels(topo_fname='./TOPO-124x124.nc',
                      wrfheight_fname='./wrfheight-124x124-22levs.nc'):
    """
    Reads the wrf height and topo files to get the asl cell height
    boundaries.
    PARAMETERS
    ----------
    topo_fname: full path to the IO/API file containing lat, lon, and
    height above sea level for STEM grid cells
    wrfheight_fname: full path to the IO/API file containing z levels
    for stem grid cells
    """
    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(topo_fname)

    f = netcdf.netcdf_file(wrfheight_fname, 'r')
    agl = f.variables['AGL'][0]  #numpy array shape(22, 124, 124): zlev, lat, lon
    f.close()

    asl = agl + topo
    return(asl)

def Haversine(p1,p2):
    """#The Haversine will calculate the distance between two coordinates.
    Input:
            p1 (list[lng,lat]) = coordinates of point 1
            p2 (list[lng,lat]) = coordinates of point 2
    Output: distance in km
    """

    # Variable setup
    lat1 = p1[0]*0.0174532925
    lat2 = p2[0]*0.0174532925
    lon1 = p1[1]*0.0174532925
    lon2 = p2[1]*0.0174532925
    R = 6371.0#km

    #d=R*math.acos(math.cos((lon1-lon2)*0.0174532925)*math.cos(lat1*0.0174532925)*math.cos(lat2*0.0174532925)+math.sin(lat1*0.0174532925)*math.sin(lat2*0.0174532925))

    dLat = (lat2-lat1)
    dLon = lon2-lon1
    a = (math.sin(dLat/2))**2 + math.cos(lat1)*math.cos(lat2)*(math.sin(dLon/2))**2
    c = 2*math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = R*c
    return d


def get_zlev(alt,y,x,zlevs):
	for level in range(22):
		print level
		print zlevs[level][y][x]
		if alt<zlevs[level][y][x]:
			return level




def main(filename):
	print filename
	# get topo data in list from bottom left to top right
	f = open('points.csv', 'r')
	topo_pts = []
	for line in f:
		temp = line.replace('\n','').split(',')
		for i in range(len(temp)):
			temp[i] = float(temp[i])
		topo_pts.append(temp)
	f.close()



	obs_pts =[]
	data = []
	extract_obs(data,obs_pts, filename)
	f = open('data2008_2009.dat','w')
	f.write(str(len(obs_pts))+'\n')
	f.write('t, x, y, z, data\n')
	zlevs = get_height_levels()


	for i in range(len(obs_pts)):
		# Calculate x and y
		d = 9999999999999999
		n=0
		for j in range(len(topo_pts)):
			distance = Haversine(topo_pts[j],[obs_pts[i][0], obs_pts[i][1]])
			if distance<d:
				n = j
				d = distance
		if n == 0:
			y = 0
		else:
			y = n/124
		x = n-y*124

		# Calculate t
		start = [0,1,1,2008]


		t = calcT(start, [obs_pts[i][2],obs_pts[i][3],obs_pts[i][4],obs_pts[i][5]])
		print t, 'dogfood'

		# Calculate p

		z = get_zlev(obs_pts[i][6],y,x,zlevs)

		# write data and print check.
		f.write(str(t)+'\t'+str(x)+'\t'+ str(y)+'\t'+str(z)+'\t'+str(data[i])+'\n')

	f.close()

def get_stem_xy(lon, lat, lon_stem, lat_stem):
    """given a set of arbitrary (lon, lat) positions, find the horizontal
    (x, y) STEM grid indices of the nearest STEM grid cell center to
    each position.
    PARAMETERS
    ----------
    lon, lat: ndarray; of arbitrary longitudes and latitudes.  Must
       contain the same number of elements.
    lon_stem, lat_stem: ndarrays; longitudes and latitudes of STEM
       grid cell centers. Must contain the same number of elements.
    """
    # convert spherical lon, lat coordinates to cartesian coords. Note
    # that these x,y,z are 3-dimensional cartesian coordinates of
    # positions on a sphere, and are thus different from the x,y,z
    # *indices* of the STEM grid.
    x, y, z = lon_lat_to_cartesian(lon, lat)
    xs, ys, zs = lon_lat_to_cartesian(lon_stem, lat_stem)

    #use a KD tree to find the nearest neighbor to (x,y,z) from the
    #points within (xs, ys, zs).  A KD tree is an data structure that
    #allows efficient queries of K-dimensional space (K here is 3).
    tree = cKDTree(np.dstack((xs.flatten(),
                              ys.flatten(),
                              zs.flatten())).squeeze())
    d, inds = tree.query(
        np.dstack((x.values,
                   y.values,
                   z.values)).squeeze(), k = 1)

    return(np.unravel_index(inds, lon_stem.shape))
    
def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R.  Taken from earthpy (TODO: note URL)
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z

if __name__ == "__main__":
    filename = os.getcwd()+'/CSV/ocs_sgp_aircraft-pfp_1_hats_event.csv'

    #main(filename)

