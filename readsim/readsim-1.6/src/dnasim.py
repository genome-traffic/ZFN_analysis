#!/usr/bin/env python

import sys
import getopt
from fractions import Fraction
import common
#import analfq


def usage(): 
  
    print >> sys.stderr, "================================================================"
    print >> sys.stderr, "         DNA simulation program v0.0"
    print >> sys.stderr, "Usage : dnasim -p 4 -h 0.05 [reference fasta] "
    print >> sys.stderr, "    --ploidy : ( default = 2 )"
    print >> sys.stderr, "    --het : heterozygosity level ( default = 0.02 that is 2%% )"
    print >> sys.stderr, "    --pre : prefix for output "


def mutate_dna( fasta, ploidy, het, pre ) :

    #common.intro( tech, ref_file, read_mu, read_dist, copy, sub_mu, in_mu, del_mu ) 
    sub_mu = het/3
    in_mu = het/3
    del_mu = het/3
    replace = "yes"

    with open(fasta) as f :
        ext = fasta.strip().split('.')[-1]
        if  len(pre) == 0 :
            prefix = ".".join( fasta.strip().split('.')[:-1] )
        else :
            prefix = pre

        for i in range(0, ploidy) :
            f.seek(0)
            fname = prefix + "." + str(i) + "." + ext
            fout = open( fname, 'w' )
            polycopy = f.readline().strip()+"_copy"+str(i)
            fout.write( polycopy + "\n")
            print >> sys.stderr, "processing " + polycopy

            if  i == 0 :
                for l in f :
                    fout.write( l )
            else :
                for l in f :
                    seq = l.strip()
                    orientation, mut_seq, mut_qual =  common.mutate( seq, "M"*len(seq), sub_mu, in_mu, del_mu, "off", replace )
                    #print >> sys.stdout, seq
                    fout.write( mut_seq + "\n" )
                    #print >> sys.stdout, mut_qual, orientation
            fout.close()


def main(): 

    ploidy = 2
    het = 0.02
    pre = ""

    try :
        if  len( sys.argv ) < 2 : 
            usage()
        else : 
            fasta = sys.argv[-1]
            opt_list, args = getopt.getopt( sys.argv[1:], 'h', ['help', 
                                                                'ploidy=', 'het=',
                                                                'pre='
                                                                 ] )
            for o, a in opt_list :
                print o, a
                if o == '--ploidy' : 
                    ploidy = int(a)
                elif  o == '--help' or o == '-h' : 
                    usage()
                elif o == '--het' : 
                    het = float(a)
                elif o == '--pre' : 
                    pre = a 
                else : 
                    print >> sys.stderr, o, a

    except getopt.error, msg:
        print >> sys.stderr,  msg


    mutate_dna( fasta, ploidy, het, pre )


###############################################################################
if  __name__ == "__main__":
    main()


