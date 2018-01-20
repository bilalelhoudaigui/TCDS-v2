import TCDS.simulation as sim
import numpy as np
import matplotlib.pyplot as plt
import dnaplotlib as dpl
# gridspec is a module which specifies the location of the subplot in the figure.
import matplotlib.gridspec as gridspec

###################### Plotting #######################

def plot_mean_sigma_genes_v2(INI_file, sigma_info, RNAPs_pos_info):

    # path to the input files (remove the "params.ini" from the path)
    path = INI_file.rpartition("/")[0] + "/"
    if path=="":
        path="."
    # read the config file
    config = sim.read_config_file(INI_file)
    # get inputs infos from the config file
    GFF_file = path+config.get('INPUTS', 'GFF')
    TSS_file = path+config.get('INPUTS', 'TSS')
    TTS_file = path+config.get('INPUTS', 'TTS')
    Prot_file = path+config.get('INPUTS', 'BARR_FIX')

    SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0')
    DELTA_X = config.getfloat('SIMULATION', 'DELTA_X')

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

    # Create the figure and all axes to draw to
    fig = plt.figure(1, figsize=(9,6)) # 3.2,2.7
    gs = gridspec.GridSpec(2, 1, height_ratios=[9, 3, 1])
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
        
    ax_mean_sig = plt.subplot(gs[0])
    ax_dna = plt.subplot(gs[1])

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
    
    # plot of sigma and mean of sigma
    plt.ion()
    ax_mean_sig.plot(cov_bp, sigma_info, linewidth= 1.5)
    ax_mean_sig.legend(loc='best', fontsize = 12)
    ax_mean_sig.set_ylim([-0.2,0.2])
    ax_mean_sig.set_xlim([0, cov_bp[-1]])

    #ax_mean_sig.set_title(r"Title goes here", fontsize = 13)
    ax_mean_sig.set_ylabel(r'Supercoiling density $(\sigma)$')
    ax_mean_sig.set_xlabel('Position (bp)')
    ax_mean_sig.plot(RNAPs_pos_info*DELTA_X, np.full(len(RNAPs_pos_info), SIGMA_0, dtype=float), 'ro', markersize=12, label = "RNA Polymerase")
    ax_mean_sig.set_ylim([-0.2, 0.2])
    plt.pause(0.001)
    plt.gcf().clear()
    plt.show()  


# Plot the supercoiling density before and after adding Gyrase (Chong experiment)
def plot_topoI_gyrase_sigma(output_dir_pre, output_dir_post):

    # get the full path
    output_dir_pre = output_dir_pre+"/all_res/save_sigma_info.npz"
    output_dir_post = output_dir_post+"/all_res/save_sigma_info.npz"
    
    # get the files
    sigma_info_pre = np.load(output_dir_pre)
    sigma_info_post = np.load(output_dir_post)
    mean_sig_WG_pre = sigma_info_pre["mean_sig_wholeGenome"]    
    mean_sig_WG_post = sigma_info_post["mean_sig_wholeGenome"]
    mean_sig_WG = np.concatenate([mean_sig_WG_pre, mean_sig_WG_post])

    # plot the result
    fig = plt.figure(1)
    plt.plot(range(0, len(mean_sig_WG)*int(2), int(2)), mean_sig_WG)
    plt.xlabel("Time (s)")
    plt.ylabel("Mean of supercoiling density")
    fig.tight_layout()
    plt.show()


# Plot the initiation rate before and after adding Gyrase (Chong experiment)
def plot_topoI_gyrase_kon(output_dir_pre, output_dir_post):

    # get the full path
    output_dir_pre = output_dir_pre+"/all_res/save_tr_info.npz"
    output_dir_post = output_dir_post+"/all_res/save_tr_info.npz"
    
    # get the files
    tr_info_pre = np.load(output_dir_pre)
    tr_info_post = np.load(output_dir_post)
    # [0,1,:] : extract from the first 3D array (0 correspond to gene0) 
    # the 2nd 2D array (1 correspond to the init_rate of the gene0)
    # and all the value during the simulation (:)  
    init_rate_pre = tr_info_pre["tr_info"][0,1,:]
    init_rate_post = tr_info_post["tr_info"][0,1,:]
    tr_info_all = np.concatenate([init_rate_pre, init_rate_post])

    # plot the result
    fig = plt.figure(1)
    plt.plot(range(0, len(tr_info_all)*int(2), int(2)), tr_info_all)
    plt.xlabel("Time (s)")
    plt.ylabel("Initiation rate")
    fig.tight_layout()
    plt.show()
