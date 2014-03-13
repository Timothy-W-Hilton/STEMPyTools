from scipy.io import netcdf
import math, numpy, os

def daysInMonth(month,leap):
	"""Returns the number of days in a month. month is an integer from 1 to 12"""
	if month==2 and leap==True:
		l = 1
	else:
		l = 0
	return {1:31, 2:28+l, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 10:31, 11:30, 12:31}[month]

#January	31
#February	28*
#March	31
#April	30
#May	31
#June	30
#July	31
#August	31
#September	30
#October	31
#November	30
#December	31

def isLeap(year):
	if year%400==0:
	   return True
	elif year%100==0:
	   return False
	elif year%4==0:
	   return True
	else:
	   return False

def calcT(start, obsT):
	"""Caclulates the stem t. Start is the start of the stem simulation
		and obsT is the observation time. Both are given as [hour, day, month, year].
	"""

	leap = isLeap(obsT[3])
	
	if obsT[3]<start[3]:
		return -1
	
	t = 0
	
	# hours contribution:
	t += obsT[0]-start[0]
	# days contribution
	t += (obsT[1]-start[1])*24
	# month contribution
	m = range(1,13,1)
	i = start[2]
	while True:
		if m[i-1] ==obsT[2]:
			break
		else:
			t+= daysInMonth(i,leap)*24
		i+=1
		if i>12:
			i=1
	# year contribution
	t+=(obsT[3]-start[3])*365*24
	return t

#def read_topo_to_csv():
#	"""Reads the lat lon of the topo file to a csv file. the point at the top is the 
#	bottom-right point and the it proceeds to the right then up.
#	"""
#	f = netcdf.netcdf_file('TOPO-124x124.nc','r')
#	variables = f.variables
#	temp1 = variables['LON']
#	temp2 = variables['LAT']

#	lon = temp1[0][0]
#	lat = temp2[0][0]


#	f2 = open('topo_coords.csv','w')
#	dim = lat.shape
#	for i in range(dim[0]):
#		for j in range(dim[1]):
#			f2.write(str(lat[i][j])+','+str(lon[i][j])+'\n')
#		
#		
#	f.close()
#	f2.close()

def get_height_levels():
	"""Reads the wrf height and topo files to get the asl cell height boundaries."""
	f = netcdf.netcdf_file('TOPO-124x124.nc','r')
	variables = f.variables
	topo = variables['TOPO'][0][0]#numpy array, shape (124, 124): lat, lon

	f.close()
	
	f= netcdf.netcdf_file('wrfheight-124x124-22levs.nc','r')
	variables = f.variables
	agl = variables['AGL'][0]  #numpy array shape(22, 124, 124): zlev, lat, lon
	# init asl variable
	asl = numpy.zeros((22,124,124))
	# assing asl values = topo + agl
	for level in range(22):
		for row in range(124):
			for col in range(124):
				asl[level][row][col] = agl[level][row][col] + topo[row][col]
	f.close()
	
	return asl
				
	


def extract_obs(data,obs_pts,filename):
	"""Extracts the observation point coords.
	lon = index 22
	lat = index 21"""

	f = open(filename,'r')
	f.readline()
	for line in f:
		s = line.split(',')
		obs_pts.append([float(s[21]),float(s[22]),int(s[4]),int(s[3]),int(s[2]),int(s[1]),float(s[23])])	
		data.append(float(s[11]))
	return 
	
	
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
			
filename = os.getcwd()+'/CSV/ocs_sgp_aircraft-pfp_1_hats_event.csv'
			
main(filename)
	
	

