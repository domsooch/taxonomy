############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################

import os
import sys
import random
import fileinput, math


    



verbose = 1

###RUN CONTROL

RunControl = {
              'labelseqs':1,
              
              }



def SaveLine(label, seq, ofp, fn):
    seqlen = len(seq)
   
    if seqlen > 1000:
        print 'seqlen: %i'%seqlen,
        x = 0
        ll = label.replace('\n', '').split('|')
        countsegment = ll[1]
        for s in range(0, seqlen, 1000):
            print " s %i"%s
            x+=1
            #ll[1] = "%i-%s"%(x, countsegment)
            label = '|'.join(ll)+"|%i-%i"%(s, s+1000)
            
            csv_line = ','.join([str(index), '', seq[s:s+1000], label, str(round(seqlen/100, 0)+1), '1',fn]) + '\n'
            #print csv_line
            ofp.writelines(csv_line)
    else:
        csv_line = ','.join([str(index), '', seq, label, str(round(seqlen/100, 0) + 1), '1',fn]) + '\n'
        ofp.writelines(csv_line)
    
if RunControl['labelseqs']:
    ud = "M:\\GGenomics\\PaulSChaud_150703\\seqdb"
    opath = os.path.join(ud, 'labelled_viruses.csv')
    ofp = open(opath, 'w')
    LL = ['index',    'seq_accession',    'seq_custom',    'targetname',
          'probes',    'replicates',    'notes']
    l = ','.join(LL)+'\n'
    ofp.writelines(l)
    fnLst = os.listdir(ud)
    index = 0
    label = 'NoLabel'
    saveLine = False
    for fn in fnLst:
        p = os.path.join(ud, fn)
        txID = fn.split('_')[0]
        seq = ''
        for line in fileinput.input(p):
            if '>' in line:
                index +=1
                line = line.replace('\n', '')
                print line
                line = line.replace(',', '')
                line = line.replace('>', '')
                lstr = line.split(' ')
                label = "%s|%s|%s\n"%(txID, lstr[0], ''.join(lstr[1:4]))
                label = label.replace('||', '|').replace('\n', '')
                print label
                if seq: saveLine = True
            else:
                seq += line.replace('\n', '')
            if saveLine:
                saveLine = False
                SaveLine(label, seq, ofp, fn)
                seq = ''
        SaveLine(label, seq, ofp, fn)
    ofp.close()
    
    
    
              
#######################
#######  CODA   #######
#######################
    
raw_input(' It appears this module ran ok:')        



