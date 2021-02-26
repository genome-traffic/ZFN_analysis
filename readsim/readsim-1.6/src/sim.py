#!/usr/bin/env python

# usage: ./cmd <fastafile> 

import sys
import random
import math
from random import shuffle
import time
import datetime

# my library
import mut
import common


MAX_LEN = 160000
MIN_LEN = 20


chr_name = ''
    
def get_chr( fasta, seq, read_len ) :

    chr_name = ''

    while True :
        line = fasta.readline().strip()
        if  line.startswith('>') :
            chr_name = line.strip().strip('>').split(' ')[0]
            #print >>  sys.stderr, chr_name
            break
        elif len(line) == 0 :
            break
        else :
            seq += line
            #print >> sys.stderr, '2', line

    return chr_name, seq

def get_more_lines( fasta, seq, read_len, base_per_line ) :

    if  len(seq) < read_len :  
        r = int( math.ceil( read_len/float(base_per_line) ) )
        for c in range(r) :
            line = fasta.readline().strip()
            seq = seq + fasta.readline().strip()
    else : 
        pass

    return seq

def readsim( sub_cmd, tech, prefix, ref_file, rev_strd, replace, read_mu, read_dist, cov_mu,  
             sub_mu, in_mu, del_mu, times ) :

    common.intro_readsim( tech, ref_file, read_mu, read_dist, cov_mu, 
                  sub_mu, in_mu, del_mu ) 

    fasta  = open(ref_file)
    orientation = ''
    i = 0
    read_num = 0

    ################################################################## 
    firstline = fasta.readline().strip()
    chr = firstline.split()[0].strip('>').split(':')[0]

    ################################################################## 
    print >> sys.stderr, "{0} {1} reading {2}".format( common.INFO_002, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), ref_file )
    seq = ''
    for line in fasta :
        seq += line.strip()

    print >> sys.stderr, "{0} {1} done.".format( common.INFO_002, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') )
    ################################################################## 
    if  read_dist == 'uniform' : 

        lens = [read_mu]

    elif  read_dist == 'normal' or read_dist == 'exp' :

        lens = []
        target = cov_mu * len( seq ) 
        s = 0

        print >> sys.stderr, "{0} {1} generating length given distribution.".format( common.INFO_002, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') )

        if  read_dist == 'normal' : 
            while  s < target :  
                l = common.get_num_from_normal( read_mu ) 
                lens.append( l ) 
                s += l

        elif  read_dist == 'exp' : 
            while  s < target :  
                l = common.get_num_from_exp( read_mu )
                lens.append( l ) 
                s += l
        else : 
            pass    
    
    else :  # with read_dist file includes read lengths

        with open(read_dist) as f:
            lens = [ int(x.strip()) for x in f.readlines() ]

        print >> sys.stderr, "{0} Randomly mixing read lengths to prevent any bias".format( common.INFO_011 )
        shuffle(lens)

    ################################################################## 
    [forward, backward] = mut.sim_reads_given_chr( sub_cmd, prefix, seq, lens, 
                                                   rev_strd, replace, read_mu, read_dist, cov_mu,
                                                   sub_mu, in_mu, del_mu, times ) 

 
    # if ends

    print >> sys.stderr, "{0} Total {1} reads are generated; {2} is forward, {3} is reversed".format( common.INFO_010, forward + backward, forward, backward )
    fasta.close()


###############################################################################
    
def qualfromfa( fasta ) : 

    fa = open( fasta, 'r' ) 
    out_fq = open( fasta.rstrip('a') + 'b', 'w' )
    out_fb = open( fasta.rstrip('a') + 'q', 'w' )

    read_name = ''
    uplmt = MAX_LEN
    lwlmt = MIN_LEN

    for line in fa : 
        
        if  line.startswith('>') : 
            read_name = line
        else :
            if  len( line.strip() ) == 0 :
                # remove zero size reads
                pass
            #elif len( line.strip() ) < lwlmt :  # since CA can't deal with reads longer than lwlmt
            #    pass
            elif len( line.strip() ) > uplmt :  # since CA can't deal with reads longer than uplmt
                out_fb.write( read_name )
                out_fb.write( line[:uplmt] + '\n' )
                out_fq.write( read_name )
                out_fq.write( 'M'*uplmt + '\n' )
            else : 
                out_fb.write( read_name )
                out_fb.write( line )
                out_fq.write( read_name )
                out_fq.write( 'M'*len( line.strip() ) + '\n' )

    out_fb.close()
    out_fq.close()
    fa.close()


###############################################################################
# muate each reads
def mutate(  sub_cmd, tech, prefix,  ref_file, replace, copy, sub_mu, in_mu, del_mu ) :

    read_mu = 0
    read_dist = ''

    common.intro( tech, ref_file, read_mu, read_dist, copy, sub_mu, in_mu, del_mu ) 

    fasta  = open(ref_file)
    seq = ''
    firstline = fasta.readline().strip()  # skip the first line
    i = 0

    if  len( firstline.split(' ') ) == 2 : 
        print >> sys.stderr, firstline
        init = int( firstline.split(' ')[1].split(':')[1].split('-')[0] ) # skip the first line '>'

    elif len( firstline.split(' ') ) == 3 : 
        init = 1
    else : 
        init = 0

    
    pos = 0
    read_num = 0

    if  sub_cmd == "fq" : 
        out_fq = open( prefix + ".fastq", 'w' )
        print >> sys.stderr, "{0} {1} are created".format( common.INFO_002, prefix + ".fastq" )
    elif  sub_cmd == "fa" : 
        out_fa = open( prefix + ".fasta", 'w' )
        print >> sys.stderr, "{0} {1} are created".format( common.INFO_003,  prefix + ".fasta" )
    #elif  sub_cmd == "fafq" :  # depreciated from v1.6
    #    out_fa = open( prefix + ".fa", 'w' )
    #    out_fq = open( prefix + ".fq", 'w' )
    #    print >> sys.stderr, "{0} {1}, {2} are created".format( common.INFO_004, prefix + ".fa", prefix + ".fq" )
    else : 
        print >> sys.stderr, "{0} sub command {1} is NOT supported".format( common.ERR_001, sub_cmd )

    total_base = 0
    i = 0

    #readlines() would be considered

    # take one long reads

    for line in fasta : 
        seq += line.strip()

    #print len(seq)
    #out_fa.write(seq + '\n')

    for i in range( 0, copy ) : 
        orientation, mut_seq, mut_qual =  common.mutate( seq, "M"*len(seq), sub_mu, in_mu, del_mu, "off", replace )
        if  sub_cmd == "fa" : 
            out_fa.write(mut_seq + '\n')
            out_fa.write(mut_qual + '\n')
        else : 
            print >> sys.stderr, "{0} sub command {1} is NOT supported".format( common.ERR_001, sub_cmd )
 

    if  sub_cmd == "fq" : 
        out_fq.close()
    elif  sub_cmd == "fa" : 
        out_fa.close()
    elif  sub_cmd == "fafq" : 
        out_fa.close()
        out_fq.close()
 
    fasta.close()


