from gridData import Grid 
import numpy as np 
import MDAnalysis as md 
#import matplotlib.pyplot as plt 


#load a gro file so that you can get the protein coordinates
prot = md.Universe('../aligned_alpha.gro')
prot_center = prot.select_atoms('protein and resid 188-230')   #used to be 205 to 238 but now we want SF and less of S6 by the C term



def opendx(x,y):
	#x = working dir, y = dx file
	g = Grid(x+y)
	return g.grid, np.asarray(g.grid.shape), g.origin, g.delta   
	## returning the grid, the shape, origin, and size of bins which will be used to scale the bins later

#this generates grid, num bins, origin of grid, and size of bins 
d = './'
num_dx = 2005
z_frame = []

for vmap in xrange(num_dx):
	f = opendx(d, str(vmap) +'.smaller.dx')#.format(vmap)

	def integrate_for_z(prot, dx):
		com = prot.center_of_mass()[0:2]
		grid = dx[0]
		origin = dx[2][0:2]
		size = dx[3][0:2]
		prot_dx = ((com - origin) / size).astype(int)  
		# if com - origin / size => this tells you the center of mass of the protein on x and y in terms of dx coordinates
		# if origin - com * size => this should be the center of the dx file in terms of the protein
		x = prot_dx[0]
		y = prot_dx[1]
		sum_z= [] 
		z_coord = []
		for z in xrange(dx[1][2]):   #this is 133 ATM, so there are 133 z points that we will search
									# the x-15 and x+15 and y-15 y+15 will give us a box of 20x 20 for this one square of z
									# so we divide by 400 to get the average values for the box
			sum_z.append(np.sum(grid[(x-15):(x+15), (y-15):(y+15), z] / 400))
			z_coord.append(z * dx[3][2] + dx[2][2])
				### getting back to coordinates by multiplying bins (z) by size of z bins (dx[3][2]) and adding the origin for reference (dx[2][2]) 
		return sum_z, z_coord 


	z_frame.append(integrate_for_z(prot_center, f))
	print vmap 


np.save('z_frame_protein_coords.npy', z_frame)

