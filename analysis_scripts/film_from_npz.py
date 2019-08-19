import sys
import numpy as np
import TSC as sim
import TSC_plotting as pl
import matplotlib
import os
matplotlib.rcParams.update({'font.size': 13})

INI_file=sys.argv[1]
output_dir=sys.argv[2]

output_dir_res = output_dir+"/all_res"

sigma_info = np.load(output_dir_res+"/save_sigma_info.npz")
RNAPs_info = np.load(output_dir_res+"/save_RNAPs_info.npz")["RNAPs_info"]
sc=pl.get_SC_array(INI_file,output_dir,compute_topoisomerase=False)

for i, rnap in enumerate(RNAPs_info[:,1,:].T):
        pl.plot_genome_and_features("test_%d"%i, INI_file, signals=[("SC","blue",sc[i])], RNAPs=rnap, width=4, height=3, hlims=(-0.2,0.2))
        os.remove("test_%d.pdf"%i)
        os.remove("test_%d.svg"%i)
os.system("ffmpeg -r 10 -i 'test_%d.png' movie.mp4")
print("output in movie.mp4")
