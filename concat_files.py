import glob
import os
import pickle 
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_folder",type=str,default=None,
                    dest="input_folder", help="name of folder in output_plots")
args = parser.parse_args()
data = []
for filename in glob.glob(args.input_folder+'/*.dat'):
    with open(os.path.join(os.getcwd(), filename), 'r') as f:
        with open(filename, 'rb') as f:
            u = pickle.Unpickler(f,encoding="latin1")
            data.append(u.load())
data = np.concatenate((data), axis =0)

pickle_out = open(args.input_folder + ".dat","wb")
pickle.dump(data,pickle_out)
pickle_out.close()
