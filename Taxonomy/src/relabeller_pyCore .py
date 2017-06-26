############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################

import os
import sys
import random
##Location
WORK = {'subRoutineDir':'c:/work/Python/New_Core/', #
        'WorkDir' : os.getcwd()+'/'}
CDrive = {'subRoutineDir':'C:\\Users\\dominic\\EclipseWorkspace\\new_pyCore\\src\\',
        'WorkDir' :os.getcwd()+'/',
        }
location = CDrive

subRoutineDir = location['subRoutineDir']

sys.path.append(subRoutineDir[:-1])





if __name__=='__main__':
    ifnLst = [
              'M:/GGenomics/20170615_design/strep_sequence.fasta', 'M:/GGenomics/20170615_design/entero_sequence.fasta' ]
    ofp = open('M:/GGenomics/20170615_design/bgdb_entero-strep.fasta', 'w')
    for ifn in ifnLst:
        ifp = open(ifn, 'r')
        buff = ifp.read(10000000)
        while buff:
            Lst = buff.split('\n')
            for i in range(len(Lst)):
                line = Lst[i]
                if '>' in line:
                    ll = line[1:].split(' ')
                    acc = ll[0]
                    group = '_'.join(ll[1:3])
                    print line, group
                    Lst[i] = '>%s|%s'%(group, acc)
            ofp.write('\n'.join(Lst))
            buff = ifp.read(10000000)
    ofp.close()

 



