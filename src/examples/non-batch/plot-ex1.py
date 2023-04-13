import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import json
import _pickle as pickle


def readStatsFile(arg_filename):
    # file open and create json object
    f1 = open(arg_filename)
    data = json.load(f1)
    f1.close()
    print(arg_filename)

    # extract site
    readings = data['readings']
    np_1d_snapshot_time = np.array([reading['time'] for reading in readings])
    np_2d_snapshot_cov = np.array([reading['species coverage'] for reading in readings])
    np_2d_snapshot_procs = np.array([reading['process counts'] for reading in readings])
    lattice = data['lattice size']
    #print(lattice)
    n_sites = lattice[0] * lattice[1] * lattice[2]


    return np_1d_snapshot_time, np_2d_snapshot_cov, np_2d_snapshot_procs, n_sites 

def data_reduce(arg_time, arg_cov, arg_procs):
    new_ints = 10000
    if (len(arg_time)-1 <= new_ints):
        return arg_time, arg_cov, arg_procs
    reduce_ratio = int((len(arg_time) - 1) / new_ints)
    print('reduce_ratio = ', reduce_ratio)

    red_time = [arg_time[i*reduce_ratio] for i in range(new_ints+1)]
    red_coverage = [arg_cov[i*reduce_ratio] for i in range(new_ints+1)]
    red_proc_counts = [arg_procs[i*reduce_ratio] for i in range(new_ints+1)]
    #print(len(red_time))
    return red_time, red_coverage, red_proc_counts

time_full, coverage_full, proc_counts_full, n_sites = readStatsFile("stats-ex1.json")
# number of species needs to be specified
n_species = len(coverage_full[0])

stats_time, stats_coverage, stats_proc_counts = data_reduce(time_full, coverage_full, proc_counts_full)


def speciesplot(arg_time, arg_cov):
    ints = len(arg_time)
    #print("ints= ", ints)
    #print('n_species= ',n_species)

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

speciesplot(stats_time, stats_coverage)

def tof_plot(arg_time, arg_procs):
    ints = len(arg_time)
    n_procs = len(arg_procs[0])
    #print(n_procs)

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
    
tof_plot(stats_time, stats_proc_counts)

def readDumpFile(arg_filename):
    # file open and create json object
    f1 = open(arg_filename)
    data = json.load(f1)
    f1.close()

    # extract coordinates
    coords = data['coordinates']

    np_2d_coords = np.array(np.array(coords['data']).reshape(coords['shape'][0], coords['shape'][1]))

    # extract site
    sites = data['sites']
    np_1d_snapshot_time = np.array([site['time'] for site in sites])
    np_2d_snapshot_sites = np.array([site['data'] for site in sites])
                            
    return np_2d_coords, np_1d_snapshot_time, np_2d_snapshot_sites


coords, dump_time, sites = readDumpFile("dump-ex1.json")
# number of species needs to be specified
n_species = 1
for i in sites:
    n_species=max(max(i),n_species)
#print('n_species= ',n_species)

def plotLattice(arg_coords):
    x = arg_coords[:,0]
    y = arg_coords[:,1]
    
    fig = plt.figure(figsize=(6,6))    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(x,y,s=100,marker = 'o')
    plt.title('Simulated Lattice Sites')
    plt.show()

plotLattice(coords)

def plotSnapshot(idx, arg_coords, arg_time, arg_sites):
    x = arg_coords[:,0]
    y = arg_coords[:,1]
    z = arg_sites[idx,:]
    s = np.array([(zval > 0)*100 for zval in z])
    #n_species=max(z)
    
    fig = plt.figure(figsize=(6,6))    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.scatter(x,y,s=s,c=z,marker = 'o', cmap = plt.cm.jet, vmin=0, vmax=n_species )
    plt.xlabel('Angstroms')
    plt.ylabel('Angstroms')
    plt.title('Frame '+ str(idx))
    plt.colorbar()
    plt.show()

# initial
#plotSnapshot(0, coords, dump_time, sites)

# second
#plotSnapshot(1, coords, dump_time, sites)

# final
plotSnapshot(len(dump_time)-1, coords, dump_time, sites)


def animateSnapshots(arg_coords, arg_time, arg_sites, arg_filename):
    x = arg_coords[:,0]
    y = arg_coords[:,1]
    z = arg_sites[0,:]
    s = np.array([(zval > 0)*100 for zval in z])
    
    
    fig = plt.figure(figsize=(12,6))
    plt.gca().set_aspect('equal', adjustable='box')
     
    sca = plt.scatter(x,y,s=s,c=z, marker = 'o', cmap = plt.cm.jet, vmin=0, vmax=n_species );
    #plt.colorbar()
    
    def update(idx):
        sca.set_sizes([(zval > 0)*100 for zval in arg_sites[idx,:]]) 
        sca.set_array(arg_sites[idx,:])
        return sca

    anim = animation.FuncAnimation(fig, update, frames=len(arg_time), interval=500)
    plt.show()

    #anim.save(arg_filename, fps=60)

animateSnapshots(coords, dump_time, sites, "ex1.mp4")