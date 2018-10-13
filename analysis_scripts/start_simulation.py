import sys
import TCDS.simulation as sim
INI_file=sys.argv[1]
try:
    output_dir=sys.argv[2]
    sim.start_transcribing(INI_file, output_dir)
except:
    sim.start_transcribing(INI_file)
