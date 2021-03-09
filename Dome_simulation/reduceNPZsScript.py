import numpy as np
import os
import sys

def main():
    
    # reduces all files in ./SIM_data
    
    for SIMname in os.listdir('./SIM_data'): # finds directories
        if os.path.isdir('./SIM_data/'+SIMname):
            fnameBase= SIMname # directory name in SIM_data

            Nsims = 0 # counts Nsims
            path = './SIM_data/'+fnameBase+'/' # path to NPZs
            for file in os.listdir(path):
                if SIMname in file: # only looks at relevant NPZs
                    if os.path.isfile(os.path.join(path, file)):
                        Nsims +=1

            for i in range(Nsims):
                with np.load(path+fnameBase+'_%03d.npz'%i) as mat: # loads in large NPZ
                    
                    # and generates smaller (excludes 4 variables)
                    fname = path+fnameBase+'_%03d.npz'%i
                    print(fname)
                    saveFunction = "np.savez(fname, "
                    for variable in mat:
                            if variable != 'borderPath' and variable != 'i_lookup' and variable != 'j_lookup' and variable != 'overlaps':
                                array = 'mat['+'"'+variable+'"'+']'
                                saveFunction = saveFunction+variable+ '='+array+','
                    saveFunction = saveFunction[:-2] + '])'
                    exec(saveFunction)   

if __name__== "__main__":
    main()