#!/usr/bin/env python

import sys
import getopt
from fractions import Fraction

import sim
import mut
#import analfq


def usage(): 
  
    print >> sys.stderr, "\nreadsim [command] [sub command] [reference fasta] <options>"
    print >> sys.stderr, "    [command] sim | mutate"
    print >> sys.stderr, "    [sub command] depends on main command"
    print >> sys.stderr, "        sim" 
    print >> sys.stderr, "            fq - generate fastq file, including reads and quality values"     
    print >> sys.stderr, "            fa - generate fasta file, including reads only"     
    print >> sys.stderr, "            fafq - (depreciated) generate fasta and fastq file separately"     
    print >> sys.stderr, "    [reference fasta] specify path and file name"
    print >> sys.stderr, "    <options> for sim"
    print >> sys.stderr, "    --pre : The prefix. The last file name will be prefix.{fast|fastq|fa|fq}"
    print >> sys.stderr, "    --rev_strd : on | off, 'on' means create backward strands as well as forward strands, randomly half and half. 'off' means no backward strands, every read is forward strand."
    print >> sys.stderr, "    --tech : The sequencing technology that you want to simulate. We do support pacbio, pacbio_ec (pacbio error corrected) and nanopore"
    print >> sys.stderr, "    --read_mu : The average on read length"
    print >> sys.stderr, "    --cov_mu : The overall coverage"
    print >> sys.stderr, "    --read_mu : The average on read length"
    print >> sys.stderr, "    --read_dist : The overall distribution of read length. Choose among {uniform, normal, exp}. You can specify file, in which each line has a length of a read. File name cannot be {uniform, normal, exp}"
    print >> sys.stderr, "    --cov_mu : The overall coverage"
    print >> sys.stderr, "    --err_sub_mu : The average on substitution rate. (ex) 0.01 means 1% of error rate"
    print >> sys.stderr, "    --err_in_mu : The average on insertion rate. (ex) 0.01 means 1% of error rate"
    print >> sys.stderr, "    --err_del_mu : The average on deletion rate. (ex) 0.01 means 1% of error rate"
    print >> sys.stderr, "    --times : used along with \"--read_dist read_lengh_file\" so that len is taken from a file. This decide whether the len is used itself or multiplied. Default is 1, which is used as it is"
    print >> sys.stderr, "    --replace : replace non[A|C|G|T] with A. Non[A|C|G|T] means any protein or N."
    print >> sys.stderr, "    <options> for mutate"
    print >> sys.stderr, "    --pre : The prefix. The final file name will be \"prefix.{fast|fastq|fa|fq}\""
    print >> sys.stderr, "    --tech : The sequencing technology that you want to simulate. We do support pacbio, pacbio_ec (pacbio error corrected) and nanopore"
    print >> sys.stderr, "    --copy : How many mutated copies?"



def set_tech(tech, read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu):
    
    read_mu_2 = 0
    read_dist_2 = ''
    cov_mu_2 = 0 
    err_sub_mu_2 = 0
    err_in_mu_2 = 0   
    err_del_mu_2 = 0

    if  tech == 'pacbio_ec' :
        if read_mu == 0       : read_mu_2 = 5000
        if read_dist == ''    : read_dist_2 = 'exp'
        if cov_mu == 0        : cov_mu_2 = 10
        if err_sub_mu == 0    : err_sub_mu_2 = 0.0033
        if err_in_mu == 0     : err_in_mu_2 = 0.0033
        if err_del_mu == 0    : err_del_mu_2 = 0.0033

    elif  tech == 'pacbio' :
        if read_mu == 0       : read_mu_2 = 5000
        if read_dist == ''    : read_dist_2 = 'exp'
        if cov_mu == 0        : cov_mu_2 = 10
        if err_sub_mu == 0    : err_sub_mu_2 = 0.01
        if err_in_mu == 0     : err_in_mu_2 = 0.12
        if err_del_mu == 0    : err_del_mu_2 = 0.02

    elif  tech == 'nanopore' :
        read_mu_2 = 100000
        read_dist_2 = 'normal'
        cov_mu_2 = 10
        err_sub_mu_2 = 0.0005
        err_in_mu_2 = 0.0005
        err_del_mu_2 = 0.0005

    elif  tech == 'perfect' :
        err_sub_mu_2 = 0
        err_in_mu_2 = 0
        err_del_mu_2 = 0
   
    elif tech == '' : 
        pass

    return (read_mu_2, read_dist_2, cov_mu_2, err_sub_mu_2, err_in_mu_2, err_del_mu_2)

 
def main(): 

    tech = ''
    rev_strd = ''
    replace = 'no'
    fastq = ''
    fasta = ''
    prefix = ''
    ref_fa = ''
    read_mu = 0
    read_dist = ''
    cov_mu = 0.0
    #cov_dist = 'poisson'
    err_sub_mu = 0.0
    #err_sub_dist = 'normal'
    err_in_mu = 0.0
    #err_in_dist = 'normal'
    err_del_mu = 0.0
    #err_del_dist = 'normal'
    times = 1
    copy = ''

    try :
        if  len( sys.argv ) <= 4 : 
            usage()

        else : 
            cmd = sys.argv[1]
            sub_cmd = sys.argv[2]
            opt_list, args = getopt.getopt( sys.argv[3:], 'h', ['ref=', 'fa=', 'pre=', 'tech=',
                                                                'rev_strd=', 
                                                                'replace=', 
                                                                'read_mu=', 'read_dist=', 
                                                                'cov_mu=', #'cov_dist=', 
                                                                'err_sub_mu=', # 'err_sub_dist=', 
                                                                'err_in_mu=', # 'err_in_dist=',
                                                                'err_del_mu=', # 'err_del_dist='
                                                                'times=',# mean*times
                                                                'copy=' # how many mutated copies
                                                               ] )
            for o, a in opt_list :
                #print o, a
                if o == '--ref' : 
                    ref_fa = a
                elif  o == '--fa' :
                    fasta = a
                elif o == '--pre' : 
                    prefix = a
                elif  o == '--tech' :
                    tech = a
                elif  o == '--copy' :
                    copy = int(a)
                elif  o == '--rev_strd' :
                    rev_strd = a
                elif  o == '--replace' :
                    replace = a
                elif o == '--read_mu' : 
                    read_mu = int(a)
                elif o == '--read_dist' : 
                    read_dist = a
                elif o == '--cov_mu' : 
                    cov_mu = float(a)
                elif o == '--err_sub_mu' :
                    if  '/' in a : 
                        err_sub_mu = float(a.split('/')[0])/float(a.split('/')[1])
                    else : 
                        err_sub_mu = float(a)
                elif o == '--err_in_mu' :
                    if  '/' in a : 
                        err_in_mu = float(a.split('/')[0])/float(a.split('/')[1])
                    else : 
                        err_in_mu = float(a)
                elif o == '--err_del_mu' :
                    if  '/' in a : 
                        err_del_mu = float(a.split('/')[0])/float(a.split('/')[1])
                    else : 
                        err_del_mu = float(a)
                elif o == '--times' :
                    times = float(a)
                else : 
                    print >> sys.stderr, "[ERROR:010] Option (", o, ") are not available"
        
            # call read simulator by tech
            #print cmd, sub_cmd
 
            if  cmd == 'sim' : 

                if  tech == 'pacbio_ec' :
                    if read_mu == 0       : read_mu = 2000
                    if read_dist == ''    : read_dist = 'exp'
                    if cov_mu == 0        : cov_mu = 10
                    #if cov_dist == ''     : cov_dist = 'poisson'
                    if err_sub_mu == 0    : err_sub_mu = 0.0033
                    #if err_sub_dist == '' : err_sub_dist = 'normal'
                    if err_in_mu == 0     : err_in_mu = 0.0033
                    #if err_in_dist == ''  : err_in_dist = 'normal'
                    if err_del_mu == 0    : err_del_mu = 0.0033
                    #if err_del_dist == '' : err_del_dist = 'normal'
        
                    sim.readsim( sub_cmd, tech, prefix, ref_fa, rev_strd, replace, read_mu, read_dist, cov_mu, 
                                 err_sub_mu, err_in_mu, err_del_mu, times  )
                elif  tech == 'pacbio' :
                    if read_mu == 0       : read_mu = 2000
                    if read_dist == ''    : read_dist = 'exp'
                    if cov_mu == 0        : cov_mu = 10
                    #if cov_dist == ''     : cov_dist = 'poisson'
                    if err_sub_mu == 0    : err_sub_mu = 0.01
                    #if err_sub_dist == '' : err_sub_dist = 'normal'
                    if err_in_mu == 0     : err_in_mu = 0.12
                    #if err_in_dist == ''  : err_in_dist = 'normal'
                    if err_del_mu == 0    : err_del_mu = 0.02
                    #if err_del_dist == '' : err_del_dist = 'normal'
        
                    sim.readsim( sub_cmd, tech, prefix, ref_fa, rev_strd, replace, read_mu, read_dist, cov_mu, 
                                 err_sub_mu, err_in_mu, err_del_mu, times )
                elif  tech == 'nanopore' :
                    if read_mu == 0       : read_mu = 10000
                    if read_dist == ''    : read_dist = 'normal'
                    if cov_mu == 0        : cov_mu = 10
                    #if cov_dist == ''     : cov_dist = 'poisson'
                    if err_sub_mu == 0    : err_sub_mu = 0.03
                    #if err_sub_dist == '' : err_sub_dist = 'normal'
                    if err_in_mu == 0     : err_in_mu = 0.03
                    #if err_in_dist == ''  : err_in_dist = 'normal'
                    if err_del_mu == 0    : err_del_mu = 0.03
                    #if err_del_dist == '' : err_del_dist = 'normal'
        
                    sim.readsim( sub_cmd, tech, prefix, ref_fa, rev_strd, replace, read_mu, read_dist, cov_mu, 
                                 err_sub_mu, err_in_mu,  err_del_mu, times )
                                 #err_sub_mu, err_sub_dist, err_in_mu, err_in_dist, err_del_mu, err_del_dist )
                elif  tech == 'perfect' :
                    #read_mu = 5000
                    #read_dist = 'uniform'
                    err_sub_mu = 0
                    err_in_mu = 0
                    err_del_mu = 0
       
                    sim.readsim( sub_cmd, tech, prefix, ref_fa, rev_strd, replace, read_mu, read_dist, cov_mu, 
                                 err_sub_mu, err_in_mu,  err_del_mu, times )
        
                elif  tech == '' : 
                    if  read_mu == 0 and read_dist == '' : 
                        print >> sys.stderr, "Must input average read length for user specific setting"
                        sys.exit()
                    if  read_dist == ''    :
                        print >> sys.stderr, "Must input reads overall distribution for user specific setting"
                        sys.exit()
                    if  cov_mu == 0        : 
                        print >> sys.stderr, "Must input average coverage for user specific setting"
                        sys.exit()
                    if  err_sub_mu == 0    : 
                        print >> sys.stderr, "Must input average substitution error rate for user specific setting"
                        sys.exit()
                    if  err_in_mu == 0     : 
                        print >> sys.stderr, "Must input average insertion error rate for user specific setting"
                        sys.exit()
                    if  err_del_mu == 0    : 
                        print >> sys.stderr, "Must input average deletion error rate for user specific setting"
                        sys.exit()
        
                    sim.readsim( sub_cmd, tech, prefix, ref_fa, rev_strd, replace, read_mu, read_dist, cov_mu, 
                                 err_sub_mu, err_in_mu, err_del_mu, times )

                else : # default tech
                    print >> sys.stderr, "{0} is not supported currently.".format( tech )
                    sys.exit()
                #analfq.readdist( fastq)
            elif  cmd == 'mutate' : 
                if  tech == 'pacbio_ec' :
                    read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu = \
                        set_tech(tech, read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu)
                    sim.mutate( sub_cmd, tech, prefix, ref_fa, replace, copy, err_sub_mu, err_in_mu,  err_del_mu )

                elif  tech == 'pacbio' :
                    print "pacbio"
                    read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu = \
                        set_tech(tech, read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu)
                    sim.mutate( sub_cmd, tech, prefix, ref_fa, replace, copy, err_sub_mu, err_in_mu,  err_del_mu )

                elif  tech == 'perfect' :
                    read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu = \
                        set_tech(tech, read_mu, read_dist, cov_mu, err_sub_mu, err_in_mu, err_del_mu)
                    sim.mutate( sub_cmd, tech, prefix, ref_fa, replace, copy, err_sub_mu, err_in_mu,  err_del_mu )

                elif  tech == '' : 
                    if  read_mu == 0 and read_dist == '' : 
                        print >> sys.stderr, "Must input average read length for user specific setting"
                        sys.exit()
                    if  read_dist == ''    :
                        print >> sys.stderr, "Must input reads overall distribution for user specific setting"
                        sys.exit()
                    if  cov_mu == 0        : 
                        print >> sys.stderr, "Must input average coverage for user specific setting"
                        sys.exit()
                    if  err_sub_mu == 0    : 
                        print >> sys.stderr, "Must input average substitution error rate for user specific setting"
                        sys.exit()
                    if  err_in_mu == 0     : 
                        print >> sys.stderr, "Must input average insertion error rate for user specific setting"
                        sys.exit()
                    if  err_del_mu == 0    : 
                        print >> sys.stderr, "Must input average deletion error rate for user specific setting"
                        sys.exit()
                    
                    gen.mutate( sub_cmd, tech, prefix, ref_fa, replace, copy, err_sub_mu, err_in_mu,  err_del_mu )
        
            else : 
                print >> sys.stderr, "Not supported command currently"

    except getopt.error, msg: 
        print >> sys.stderr,  msg


###############################################################################
if  __name__ == "__main__":
    main()


