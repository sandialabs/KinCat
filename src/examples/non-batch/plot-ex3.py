import numpy as np 
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def read_dump(arg_file):
	with h5py.File(arg_file, 'r') as f:
		coord_data = f['Dump']['coordinates.data']
		x=[]
		y=[]
		n_sites = int(len(coord_data)/2)
		for i in range(n_sites):
			x.append(coord_data[i*2])
			y.append(coord_data[(i*2)+1])
		n_snapshots = f['Dump']['number of snapshots'][()]
		sites_data = []

		for n in range(n_snapshots):
			occ_data = f['Dump']['snapshot.'+str(n)+'.data'][:]
			temp_data = []
			for i in range(n_sites):
				temp_data.append(occ_data[i])
			sites_data.append(np.array(temp_data))
			
		return x,y,sites_data

def plot_lattice(x,y):
	with h5py.File(filename, 'r') as f:
		
		plt.plot(x,y,'o')
		plt.title('Lattice')
		plt.xlabel('Angstroms')
		plt.ylabel('Angstroms')
		plt.show()

def plot_snapshot(n, x, y):
	with h5py.File(filename, 'r') as f:
		occ_data = f['Dump']['snapshot.'+str(n)+'.data'][()]
		s = np.array([(zval > 0)*100 for zval in occ_data])

		plt.scatter(x,y, s=s, c=occ_data, marker = 'o', cmap=plt.cm.jet, vmin=1, vmax=n_species-1 )
		plt.colorbar()
		plt.title('Snapshot '+str(n))
		plt.xlabel('Angstroms')
		plt.ylabel('Angstroms')
		plt.show()

def animateSnapshots(x,y, sites_data):
	with h5py.File(filename, 'r') as f:
		n_snapshots = f['Dump']['number of snapshots'][()]
		print('n_snapshots=', n_snapshots)
		s = np.array([(zval > 0)*100 for zval in sites_data[0]])

		fig = plt.figure(figsize=(12,6))
		plt.gca().set_aspect('equal', adjustable='box')

		sca = plt.scatter(x,y,s=s,c=sites_data[0], marker = 'o', cmap = plt.cm.jet, vmin = 1, vmax=n_species-1)
		plt.colorbar()

		def update(idx):
			sca.set_sizes([(zval > 0)*100 for zval in sites_data[idx]])
			sca.set_array(sites_data[idx])
			return sca

		anim = animation.FuncAnimation(fig,update, frames=len(sites_data), interval = 500)
		plt.show()

filename = 'dump-ex3.hdf5'
[x_coords, y_coords, sites_data] = read_dump(filename)
n_species = 3
#for i in range(len(sites_data)):
#	n_species = np.max(np.max(sites_data[i]),n_species)

animateSnapshots(x_coords,y_coords, sites_data)

### Stats file read

def readStatsFileHDF5(arg_file):
	with h5py.File(arg_file, 'r') as f:
		# extract site
		n_snapshots_stats = np.array(f['Stats']['number of snapshots'])
		n_basis = np.array(f['Stats']['lattice.shape.basis'])
		dim1 = np.array(f['Stats']['lattice.shape.dim1'])
		dim2 = np.array(f['Stats']['lattice.shape.dim2'])
		lattice = [dim1, dim2, n_basis]
		n_species = np.array(f['Stats']['snapshot.0.coverage.Length'])-1
		n_procs = np.array(f['Stats']['snapshot.0.processcounts.Length'])

		np_1d_snapshot_time = [0 for _ in range(n_snapshots_stats)]
		np_2d_snapshot_cov = []#np.array([[0 for i in range(n_species+1)] for _ in range(n_snapshots_stats)])
		np_2d_snapshot_procs = np.array([[0 for i in range(n_procs)] for _ in range(n_snapshots_stats)])
		np_3d_snapshot_sitecov = []#np.array([[[0 for j in range(n_species+1)] for i in range(n_basis)] for _ in range(n_snapshots_stats)])
		n_sites = dim1*dim2*n_basis;

		for i in range(n_snapshots_stats):
			temp_str = 'snapshot.' + str(i)
			np_1d_snapshot_time[i] = np.copy(np.array(f['Stats'][temp_str+'.time']))
			temp_array = np.copy(np.array(f['Stats'][temp_str+'.coverage']))
			np_2d_snapshot_cov.append(np.copy(temp_array))
			np_2d_snapshot_procs[i] = f['Stats'][temp_str+'.processcounts']
			temp_array = f['Stats'][temp_str+'.speciescoverage']
			np_3d_snapshot_sitecov.append(np.copy(temp_array))
			#for j in range(n_basis):
			#	np_3d_snapshot_sitecov[i][j] = temp_array[(j*(n_species+1)):((j+1)+n_species+1)-1]

	return np_1d_snapshot_time, np_2d_snapshot_cov, np_2d_snapshot_procs, n_sites, np_3d_snapshot_sitecov 


filename_stats = 'stats-ex3.hdf5'
time_full, coverage_full, proc_counts_full, n_sites, sitecov_full = readStatsFileHDF5(filename_stats)

def speciesplot(arg_time, arg_cov):
    ints = len(arg_time)

    y0=[]
    y1=[]
    y2=[]
    for j in range(ints):
        y0.append(arg_cov[j][0])
        y1.append(arg_cov[j][1])
        y2.append(arg_cov[j][2])

    plt.plot(arg_time,y0, label="Vacant")

    plt.plot(arg_time,y1, label="O")

    plt.plot(arg_time,y2, label="CO")
    plt.legend()
    plt.title('Species Coverage')
    plt.show()

speciesplot(time_full, coverage_full)

def tof_plot(arg_time, arg_procs):
    ints = len(arg_time)
    n_procs = len(arg_procs[0])

    tofs = [[] for _ in range(ints)]
    tofs[0] = 0
    for i in range(1,ints):
        tof = 0
        for j in [18,19,20,21]:
            tof += arg_procs[i][j] - arg_procs[i-1][j]
        t_step = arg_time[i]-arg_time[i-1]
        tofs[i] = tof/(t_step*n_sites)

    plt.plot(arg_time,tofs)
    plt.title('Turn-Over Frequency for CO Oxidation')
    plt.show()
    
tof_plot(time_full, proc_counts_full)

