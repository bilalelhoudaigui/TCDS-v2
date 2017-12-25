import sys
import numpy as np
import TCDS.simulation as sim
import TCDS.plotting as plotting
import matplotlib
matplotlib.rcParams.update({'font.size': 13})

INI_file=sys.argv[1]
output_dir=sys.argv[2]

GFF_file, TSS_file, TTS_file, Prot_file, m, sigma_t, epsilon, SIGMA_0, DELTA_X, DELTA_T, RNAPS_NB, ITERATIONS_NB, OUTPUT_STEP, GYRASE_CONC, TOPO_CONC, TOPO_CTE, GYRASE_CTE, TOPO_EFFICIENCY, k_GYRASE, x0_GYRASE, k_TOPO, x0_TOPO = sim.read_config_file_v2(INI_file)
tss = sim.load_tab_file(TSS_file)
Kon = tss['TSS_strength'].values

output_dir_res = output_dir+"/withSC_Kon_%.06f/RNAPn_%s/Sig0_%s/Gyrase_%s_TopoI_%s/all_res" %(Kon[0], RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC)

sigma_info = np.load(output_dir_res+"/save_sigma_info.npz")
RNAPs_info = np.load(output_dir_res+"/save_RNAPs_info.npz")

Barr_sigma_info = sigma_info["Barr_sigma_info"]
Dom_size_info = sigma_info["Dom_size_info"]

for i, Barr_sigma_val in enumerate(sigma_info["Barr_sigma_info"]):
	one_sigma_info = np.repeat(Barr_sigma_val, sigma_info["Dom_size_info"][i])
	RNAPs_pos_info = RNAPs_info["RNAPs_info"][:, 1, i]
	plotting.plot_mean_sigma_genes_v2(INI_file, one_sigma_info, RNAPs_pos_info)
