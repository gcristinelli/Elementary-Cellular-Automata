import numpy as np
import time
import random as rd
import matplotlib.pyplot as plt

start= time.time()

np.random.seed(123)

all_lmd=[1]
#N2=np.arange(0.4,0.55,0.05)
N3=np.arange(0.6,0.75,0.05)
N4=np.arange(0.75,0.84,0.01)
N5=np.arange(0.84,1.01,0.02)

#all_lmd=np.concatenate((N1,N3,N4,N5),axis=None).tolist() #list of all the lambdas that we will consider

T=2000 #number of iteration
N=2000#number of cells
n_runs=10 #number of iteration on the same  lambda
rule1=0 #with probability (1-lmd)
rule2=184 # second rule number with probability (lmd)
rule1_bin = [int(x) for x in np.binary_repr(rule1, width=8)] # binary representation of f2 rule
rule2_bin = [int(x) for x in np.binary_repr(rule2, width=8)] # binary representation of f2 rule
rule_set = np.array([[1,1,1],[1,1,0],[1,0,1],[1,0,0],[0,1,1],[0,1,0],[0,0,1],[0,0,0]])

def ruleFunction(x,rule_bin,rule_set):
	Outx = np.zeros(np.shape(x))
	if np.sum(rule_bin)!=0:			
		# input a np array of size(n,1)
		neighArray = np.zeros((len(x),3))
		neighArray[:,0] =  np.roll(x,-1) #neighbour to the left
		neighArray[:,1] = x
		neighArray[:,2] = np.roll(x,1) #neighbour to the right
		bins = np.where(rule_bin)
		for i in range(np.sum(rule_bin)):
			Outx = Outx + (neighArray[:] == rule_set[bins[0][i]]).all(axis=1)
	return Outx
	
def density(lmd,n):
	grid = np.zeros((T,n)) # T x n grid (matrix) of zeroions (white cells)
	grid[0, :]=(np.random.rand(n)<0.3)# random first row with 50% chance of either 0 or 1
	ArrayRand = (np.random.rand(T-1,n)<lmd)
	dens=np.zeros(T)
	dens[0]=np.sum(grid[0,:])/n
	for t in range(1,T):
		grid[t] = np.abs(ArrayRand[t-1]-1)*ruleFunction(grid[t-1],rule1_bin,rule_set) + ArrayRand[t-1]*ruleFunction(grid[t-1],rule2_bin,rule_set)
		dens[t]=np.sum(grid[t,:])/n
	if lmd==1:
		plt.imshow(grid, cmap='hot', interpolation='nearest')
		plt.show()
	return dens
	
def plot_density(d_avg, d_std):
	plt.figure(3,figsize=[10,10])
	plt.grid(b=True)
	ax = plt.gca()
	ax.set_ylim([-0.01, 1])
	#plt.title("density plot for rules %d and %d with average from %d runs" % (rule1,rule2,n_runs),fontsize=18)
	plt.xlabel("$\\lambda$",fontsize=25)
	plt.ylabel("density",fontsize=25)	
	plt.errorbar(all_lmd,d_avg,yerr=d_std,fmt='.',ls=':', color='blue',ecolor='green', elinewidth=1.5, capsize=2)
	plt.savefig("density_plot_rule=%d_n=%d_t=%d_nruns=%d.png" % (rule2,N,T,n_runs), delimiter=",")
	
	
def plot_time_density (d_avg, d_std):
    plt.figure(1,figsize=[10,10])
    plt.grid(b=True)
    #plt.title("Time-density for $\\lambda$-s for rules %d and %d\n Average from %d runs" % (rule1,rule2,n_runs),fontsize=18)
    plt.ylabel("$density$",fontsize=25)  
    plt.xlabel("time t",fontsize=25)
    ax = plt.gca()
    ax.set_ylim([0, 1])
    T_plot = np.linspace(0,T,T)
    for i in range(0,len(all_lmd)):
        plt.errorbar(T_plot[::3],d_avg[i,::3],yerr=d_std[i,::3],marker='.',ls=':',elinewidth=0.5,capsize=0)
    ax.legend(all_lmd,loc='upper center', bbox_to_anchor=(0.5, 1.1),ncol=4, fancybox=True, shadow=True,fontsize=20)
    plt.savefig("density_time_rule=%d_n=%d_t=%d_nruns=%d.png" % (rule2,N,T,n_runs), delimiter=",")
 
    
def plot_volume_density ():
    d=np.zeros((n_runs,len(all_lmd),N))
    for i in range(100,N,100):
    	print("Evaluating for %d cells" % (i))
    	for j in np.arange(0,len(all_lmd)):
    		all_lmd[j] = round(all_lmd[j], 4) #round the lambdas so that it is simpler to plot 
    		print("  Evaluating %d/%d lambdas" % (j+1,len(all_lmd)))
    		for k in np.arange(0,n_runs):
    			print("   Evaluating %d/%d runs" % (k+1,n_runs))
    			d[k,j,i]=density(all_lmd[j],i)[T-1]
    d_avg=np.average(d,axis=0)
    d_std=np.std(d,axis=0)
    plt.figure(1,figsize=[16,9])
    plt.grid(b=True)
    plt.title("volume-density for $\\lambda$-s for rule $rule$=%d\n Average from %d runs" % (rule2,n_runs),fontsize=18)
    plt.ylabel("$density$",fontsize=25)  
    plt.xlabel("number of cells",fontsize=25)
    T_plot = np.linspace(0,N,N)
    for i in range(0, len(all_lmd)):
        plt.errorbar(T_plot[::100],d_avg[i,::100],yerr=d_std[i,::100],marker='.',ls=':',elinewidth=0.5,capsize=0)
    plt.legend(all_lmd,ncol=8,fontsize=18)
    plt.savefig("density_volume_rule=%d_t=%d_nruns=%d.png" % (rule2,T,n_runs), delimiter=",",dpi=250)
    
    	
    	
d=np.zeros((n_runs,len(all_lmd),T))
for j in np.arange(0,len(all_lmd)):
    all_lmd[j] = round(all_lmd[j], 4) #round the lambdas so that it is simpler to plot 
    print("Evaluating %d/%d lambdas" % (j+1,len(all_lmd)))
    for i in np.arange(0,n_runs):
    	print("     Evaluating %d/%d runs" % (i+1,n_runs))
    	d[i,j,:]=density(all_lmd[j],N) 
    	
d_avg=np.average(d,axis=0)
d_std=np.std(d,axis=0)

res=np.vstack((all_lmd,d_avg[:,T-1],d_std[:,T-1]))

np.savetxt("rule=%d_n=%d_t=%d_nruns=%d_results.csv" % (rule2,N,T,n_runs), res, delimiter=",")

plot_density(d_avg[:,T-1],d_std[:,T-1]) 
plot_time_density(d_avg, d_std)

#plot_volume_density()


	
			
                



end=time.time()

print(round((end-start)/60,3), "minutes")







	
