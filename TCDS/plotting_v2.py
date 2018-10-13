import simulation as sim
import numpy as np
import matplotlib.pyplot as plt
import dnaplotlib as dpl
# gridspec is a module which specifies the location of the subplot in the figure.
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib.patches as pat
import math
import matplotlib.font_manager as font_manager
from scipy.optimize import fsolve
plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'legend.fontsize': 9})
#plt.rcParams.update({'mathtex.fontset': "cm"})
plt.rcParams.update({'font.family': "Arial"})

global exts
exts=[".pdf",".svg",".png"]



###################### Plotting #######################




def get_cov_bp(INI_file):
    """
    Analyzes initiation file to get the array of positions in reduced units, for plotting
    """
    path = INI_file.rpartition("/")[0]
    if path=="":
        path="."
    path+="/"
    # read the config file
    config = sim.read_config_file(INI_file)
    # get inputs infos from the config file
    GFF_file = path+config.get('INPUTS', 'GFF')
    DELTA_X = config.getfloat('GLOBAL', 'DELTA_X')
    # To draw the beautiful genes we need to read the GFF, TSS and TTS files to get some info ;)
    gff_df_raw = sim.load_gff(GFF_file)
    # to get the cov_bp (a verifier)
    genome_size = sim.get_genome_size(gff_df_raw)
    genome = int(genome_size/DELTA_X)
    #print(genome)
    cov_bp = np.arange(0, genome_size, DELTA_X)
    #print(cov_bp,len(cov_bp))
    #cov_bp = np.resize(cov_bp, genome)
    return cov_bp

def plot_genome(ax_dna, INI_file):
    """
    General Function that plots a genome from an INI file and puts it into a subplot
    """
    
    # path to the input files (remove the "params.ini" from the path)
    path = INI_file.rpartition("/")[0]
    if path=="":
        path="."
    path+="/"
    # read the config file
    config = sim.read_config_file(INI_file)
    # get inputs infos from the config file
    GFF_file = path+config.get('INPUTS', 'GFF')
    TSS_file = path+config.get('INPUTS', 'TSS')
    TTS_file = path+config.get('INPUTS', 'TTS')
    Prot_file = path+config.get('INPUTS', 'BARR_FIX')

    SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    DELTA_X = config.getfloat('GLOBAL', 'DELTA_X')

    # load and get BARR_FIX positions
    prot = sim.load_tab_file(Prot_file)
    BARR_FIX = (prot['prot_pos'].values).astype(int)

    # To draw the beautiful genes we need to read the GFF, TSS and TTS files to get some info ;)
    gff_df_raw = sim.load_gff(GFF_file)
    # to get the cov_bp (a verifier)
    genome_size = sim.get_genome_size(gff_df_raw)
    genome = int(genome_size/DELTA_X)
    cov_bp = np.arange(0, genome_size, DELTA_X)
    cov_bp = np.resize(cov_bp, genome)

    gff_df = sim.rename_gff_cols(gff_df_raw)

    tss = sim.load_tab_file(TSS_file)
    Kon = tss['TSS_strength'].values

    tts = sim.load_tab_file(TTS_file)
    Poff = tts['TTS_proba_off'].values

    strands = sim.str2num(gff_df['strand'].values)

    # Color maps for formatting
    col_map = {}
    col_map['red']     = (0.95, 0.30, 0.25)
    col_map['green']   = (0.38, 0.82, 0.32)
    col_map['blue']    = (0.38, 0.65, 0.87)
    col_map['orange']  = (1.00, 0.75, 0.17)
    col_map['purple']  = (0.55, 0.35, 0.64)
    col_map['yellow']  = (0.98, 0.97, 0.35)
    col_map['grey']    = (0.70, 0.70, 0.70)
    col_map['dark_grey'] = (0.60, 0.60, 0.60)
    col_map['light_grey'] = (0.9, 0.9, 0.9)

    # CDS formatting options
    opt_CDSs = []

    Ps = []
    CDSs = []
    Ts = []

    design = []

    for i in gff_df.index.values:
        opt_CDSs.append({'label':'Gene%s \n%.03f'%(str(i+1),Kon[i]), 
                         'label_style':'italic', 
                         'label_y_offset':-5,
                         'label_size':9,
                         'color':col_map['orange']})
        # Design of the construct
        if strands[i] == True:
            # Promoters
            Ps.append({'type':'Promoter', 'name':'P%s'%str(i+1), 'start':tss['TSS_pos'][i], 
                       'end':tss['TSS_pos'][i]+5, 'fwd':strands[i], 'opts':{'color':col_map['green']}})
            # Coding Sequence
            CDSs.append({'type':'CDS', 'name':'CDS%s'%str(i+1), 'start':gff_df['start'][i],  
                         'end':gff_df['end'][i], 'fwd':gff_df['strand'][i], 'opts':opt_CDSs[i]}) 
        else:
            # Promoters
            Ps.append({'type':'Promoter', 'name':'P%s'%str(i+1), 'start':tss['TSS_pos'][i], 
                       'end':tss['TSS_pos'][i]-5, 'fwd':strands[i], 'opts':{'color':col_map['green']}})
            # Coding Sequence
            CDSs.append({'type':'CDS', 'name':'CDS%s'%str(i+1), 'start':gff_df['end'][i],  
                         'end':gff_df['start'][i], 'fwd':gff_df['strand'][i], 'opts':opt_CDSs[i]}) 
        # Terminators
        Ts.append({'type':'Terminator', 'name':'T%s'%str(i+1), 'start':tts['TTS_pos'][i], 
              'end':tts['TTS_pos'][i]+5, 'fwd':strands[i], 'opts':{'color':col_map['red']}})

        # A design is merely a list of parts and their properties
        if strands[i] == True:
            design.append(Ps[i])
            design.append(CDSs[i]) 
            design.append(Ts[i])
        else:
            design.append(Ts[i])
            design.append(CDSs[i]) 
            design.append(Ps[i])

    # Redender the DNA
    dr = dpl.DNARenderer(scale=7, linewidth=1)
    start, end = dr.renderDNA(ax_dna, design, dr.trace_part_renderers())

    # Set bounds and display options for the DNA axis
    dna_len = end-start
    ax_dna.set_xlim([cov_bp[0], cov_bp[-1]]) #start-50
    ax_dna.set_ylim([-8,8])
    #ax_dna.plot(5000, 'ro', markersize=15)
    for xc in BARR_FIX:
        ax_dna.axvline(x=xc, ymin=0.40, ymax=0.60, color='k', linewidth=5)
    
    ax_dna.plot([cov_bp[0],cov_bp[-1]], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
    ax_dna.axis('off')
    return SIGMA_0, DELTA_X, BARR_FIX, cov_bp




# --------------------------------
# Plotting functions for one timepoint, one genome


def compute_superc_distrib(Barr_pos, SC, cov_bp):
    """
    Utility function
    Computes the array of SC from Barrier positions and SC values, and arange of genome length
    """
    #print(Barr_pos)
    #print(SC)
    n=len(cov_bp)
    if len(Barr_pos)>1:
        sizes=[Barr_pos[0]]+list(Barr_pos[1:]-Barr_pos[:(-1)])+[n-Barr_pos[-1]]
        SC=[SC[-1]]+list(SC)
        #print(Barr_pos,n)
        #print(sizes)
        #print(SC)np.repeat(SC, sizes)
        #print("le",len(np.repeat(SC, sizes)))
        return np.repeat(SC, sizes)
    elif len(SC)==1:
        return np.ones(n)*SC[0]
    else:
        print("problem: inconsistent barriers and SC values...")
        print(SC,Barr_pos)
        return 1

    
def plot_superc_distrib(ax, Barr_pos, SC, cov_bp, DELTA_X, Barr_fix):
    """
    Utility function
    Computes the array of SC from Barrier positions and SC values, and arange of genome length
    """
    n=len(cov_bp)
    for i,b in enumerate(Barr_pos[:(-1)]):
        x=np.arange(b,Barr_pos[i+1]+1)*DELTA_X
        ax.plot(x,np.full(len(x),SC[i]),color="blue")
    x=np.arange(Barr_pos[0]+1)*DELTA_X
    ax.plot(x,np.full(len(x),SC[-1]),color="blue")
    x=np.arange(Barr_pos[-1],n)*DELTA_X
    ax.plot(x,np.full(len(x),SC[-1]),color="blue")
    #for b in Barr_pos:
        #ax.axvline(b*DELTA_X,color="black")
        

    
def plot_genome_and_features(outfile, INI_file, signals=None, RNAPs=None):
    """
    Plots a genome into an output figure file. Optionally, make a second plot with one or several signals along the genome and/or RNAP positions. 
    - the signals must have the same size as the genome in reduced units. They are a list of tuples (label, array of values) of genome size OR tuple (label, value, Barr_pos) to draw the distribution at given timepoint from simul output
    - the RNAPs are shown as red circles. it is an array of positions
    """
    # Create the figure and all axes to draw to

    if signals==None:
        fig = plt.figure(figsize=(7,1.5)) # 3.2,2.7
        ax_dna=plt.subplot()
        SIGMA_0, DELTA_X, BARR_FIX, cov_bp=plot_genome(ax_dna, INI_file)
        if RNAPs!=None:
            ax_dna.plot(RNAPs*DELTA_X, np.full(len(RNAPs), 0.5, dtype=float), 'o', markersize=10, color="blue", zorder=100)

    else:
        fig = plt.figure(figsize=(7,3)) # 3.2,2.7
        gs = gridspec.GridSpec(2, 1, height_ratios=[1,1])
        
        ax_sig = plt.subplot(gs[1])
        ax_dna = plt.subplot(gs[0])
        
        SIGMA_0, DELTA_X, BARR_FIX, cov_bp=plot_genome(ax_dna, INI_file)
        
        
        # plot of signals
        if signals!=None:
            for si in signals:
                if len(si)==2:
                    # case where we plot an array of values along the genome
                    slab,s=si
                    ax_sig.plot(cov_bp,s,linewidth= 1.5, label=slab)
                elif len(si)==3:
                    # case where we compute the SC distribution along the genome at a timepoint
                    slab,SC,Barr_pos=si
                    plot_superc_distrib(ax_sig,Barr_pos,SC,cov_bp,DELTA_X,BARR_FIX)
                ax_sig.legend(loc='best', fontsize = 12)
                #ax_sig.set_ylim([-0.2,0.2])
            ax_sig.set_xlim([0, cov_bp[-1]])
            ax_sig.set_ylabel(r'$\sigma(x)$')
            ax_sig.set_xlabel('Position (bp)')
            if RNAPs!=None:
                ax_dna.plot(RNAPs*DELTA_X, np.full(len(RNAPs), 0, dtype=float), 'o', markersize=10, color="blue", zorder=100)
                for x in RNAPs:
                    #print(x*DELTA_X)
                    ax_dna.axvline(x=x*DELTA_X,ymin=-1.5,ymax=0.5,color="blue",ls="--",lw=.8,zorder=110,clip_on=False)
                    #ax_dna.plot([x*DELTA_X,x*DELTA_X],[0.5,0],zorder=120)
                    # con=pat.ConnectionPatch(xyA=(x*DELTA_X, 0.5), xyB=(x*DELTA_X, 0.), coordsA="data", coordsB="data", axesA=ax_dna, axesB=ax_sig, color="red")
                    # ax_sig.add_artist(con)
            for x in BARR_FIX:
                ax_dna.axvline(x=x,ymin=-1.5,ymax=0.5,color="black",ls="--",lw=0.8,zorder=110,clip_on=False)
    plt.tight_layout()
    for ext in exts:
        plt.savefig(outfile+ext)




        
# --------------
# Analysis functions that generate tables of SC and/or k_on values from output dir

def get_SCprofiles_from_dir(output_dir,compute_topoisomerase=False,timepoints=None):
    """
    Provides a list with successive tuples of Barr_fix,SC_profile that can be used to draw the distribution. 
    - if compute_topoisomerase, then also lists for those: then the argument must be the input file: params.ini
    - timepoints is an array of indexes, from 0 to the maximal timeindex
    """
    sigma_info = np.load(output_dir+"/all_res/save_sigma_info.npz")
    RNAPs_info = np.load(output_dir+"/all_res/save_RNAPs_info.npz")
    Barr_pos = sigma_info["save_Barr_pos"]
    dom_sigma_info = sigma_info["dom_sigma_info"]
    #print(dom_sigma_info[96:110])
    # select timepoints
    if timepoints is None:
        timepoints=np.arange(len(dom_sigma_info))
        sigma=dom_sigma_info
        barr=Barr_pos
        RNAPs_pos_info = RNAPs_info["RNAPs_info"][:, 1, :]
    else:
        sigma=dom_sigma_info[timepoints]
        barr=Barr_pos[timepoints]
        RNAPs_pos_info = RNAPs_info["RNAPs_info"][:, 1, timepoints]    
    # compute topoisomerases?
    if not compute_topoisomerase:
        return [(barr[i],s) for i,s in enumerate(sigma)]
    else:
        inf=compute_topoisomerase
        config = sim.read_config_file(inf)
        # get promoter values from the config file
        m = config.getfloat('GLOBAL', 'm')
        sigma_t = config.getfloat('GLOBAL', 'sigma_t')
        epsilon = config.getfloat('GLOBAL', 'epsilon')
        # get topoisomerase concentrations
        GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC')
        TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC')
        # topoisomerase behavior
        TOPO_CTE = config.getfloat('TOPOISOMERASES', 'TOPO_CTE')
        GYRASE_CTE = config.getfloat('TOPOISOMERASES', 'GYRASE_CTE')
        k_GYRASE = config.getfloat('TOPOISOMERASES', 'k_GYRASE')
        x0_GYRASE = config.getfloat('TOPOISOMERASES', 'x0_GYRASE')
        k_TOPO = config.getfloat('TOPOISOMERASES', 'k_TOPO')
        x0_TOPO = config.getfloat('TOPOISOMERASES', 'x0_TOPO')
        # compute topo activity in 
        gyr_act=[GYRASE_CONC*1/(1+np.exp(-k_GYRASE*(s-x0_GYRASE)))*GYRASE_CTE for s in sigma]
        topo_act=[TOPO_CONC*1/(1+np.exp(k_TOPO*(s-x0_TOPO)))*TOPO_CTE for s in sigma]
        return [(barr[i],s) for i,s in enumerate(sigma)], gyr_act, topo_act


def get_SC_array(init_file, output_dir,compute_topoisomerase=False,timepoints=None):
    # same as last function except that output is a Numpy array with values at each position rather than a list of domain
    # this is helpful if you want to draw the distribution of SC or topo activity along the genome
    """
    Analyzes an initiation file and associated output directory, and computes a NumPy array with values of SC at all positions and all timepoints (matrix of size timepoints x cov_bp). 
    Arguments: 
    - Initiation file for genome description
    - Output dir for simulation data
    Options:
    - compute_topoisomerase = True: in addition to the SC matrix, computes the topoisomerase activities at all positions and timepoints (by taking the topoisomerase parameters from the initfile and applying them on the latter matrix). 
    - timepoints (= NumPy array of timepoints): restrict on a subset of timepoints
    Output: 
    NumPy array or tuple of NumPy arrays of size (genome size * time)
    """
    cov_bp=get_cov_bp(init_file)
    if not compute_topoisomerase:
        bs=get_SCprofiles_from_dir(output_dir,compute_topoisomerase,timepoints)
        return np.array([compute_superc_distrib(bsi[0], bsi[1], cov_bp) for bsi in bs])
    else:
        bs,gy,to=get_SCprofiles_from_dir(output_dir,compute_topoisomerase=init_file,timepoints=timepoints)
        #print(gy)
        sc=np.array([compute_superc_distrib(bsi[0], bsi[1], cov_bp) for bsi in bs])
        gyr=np.array([compute_superc_distrib(bs[i][0], g, cov_bp) for i,g in enumerate(gy)])
        topo=np.array([compute_superc_distrib(bs[i][0], t, cov_bp) for i,t in enumerate(to)])
        return sc,gyr,topo



def SC_numerical_solution(GYRASE_CONC,TOPO_CONC,GYRASE_CTE=0.01,TOPO_CTE=0.005,k_GYRASE=50,k_TOPO=80,x0_GYRASE=.016,x0_TOPO=-.04):
    """
    Computes the equilibrium SC value from topoisomerase concentrations. 
    """
    func = lambda sig0 : -GYRASE_CONC*1/(1+np.exp(-k_GYRASE*(sig0-x0_GYRASE)))*GYRASE_CTE + TOPO_CONC*1/(1+np.exp(k_TOPO*(sig0-x0_TOPO)))*TOPO_CTE
    sig0_initial_guess = -0.03
    sig0 = fsolve(func, sig0_initial_guess)[0]
    return sig0


def plot_promoter_response_and_SCvalues(INI_file,outfile=None):
    """
    For given simulation, plot the promoter response curve together with initiatil and equilibrium SC values
    """
    if outfile is None:
        outfile=INI_file.split(".")[0]+"_promoter"
    config = sim.read_config_file(INI_file)
    # get promoter values from the config file
    m = config.getfloat('GLOBAL', 'm')
    sigma_t = config.getfloat('GLOBAL', 'sigma_t')
    epsilon = config.getfloat('GLOBAL', 'epsilon')
    # get topoisomerase constants
    GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC')
    TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC')
    TOPO_CTE = config.getfloat('TOPOISOMERASES', 'TOPO_CTE')
    GYRASE_CTE = config.getfloat('TOPOISOMERASES', 'GYRASE_CTE')
    # topoisomerase behavior
    k_GYRASE = config.getfloat('TOPOISOMERASES', 'k_GYRASE')
    x0_GYRASE = config.getfloat('TOPOISOMERASES', 'x0_GYRASE')
    k_TOPO = config.getfloat('TOPOISOMERASES', 'k_TOPO')
    x0_TOPO = config.getfloat('TOPOISOMERASES', 'x0_TOPO')
    # sigma0,sigma_eq
    sigma_eq=SC_numerical_solution(GYRASE_CONC,TOPO_CONC,GYRASE_CTE,TOPO_CTE,k_GYRASE,k_TOPO,x0_GYRASE,x0_TOPO)
    try:
        SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    except:
        SIGMA_0=sigma_eq
    # -------------------------
    prom = lambda sig: np.exp((1/(1+np.exp((sig-sigma_t)/epsilon)))*m)
    #
    fig = plt.figure(figsize=(4,3)) # 3.2,2.7
    sigs=np.arange(-.12,.04,.005)
    plt.plot(sigs,prom(sigs),color="black")
    plt.axvline(SIGMA_0,color="gray",ls="--",lw=.5,label="initial")
    plt.axvline(sigma_eq,color="gray",lw=.5,label="equil")
    plt.xlabel("σ")
    plt.ylabel("supercoiling activation factor")
    plt.legend()
    plt.tight_layout()
    for ext in exts:
        plt.savefig(outfile+ext)
    plt.close()
    
