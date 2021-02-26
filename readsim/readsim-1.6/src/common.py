#!/usr/bin/env python

# usage: ./cmd <fastafile> 

import sys
import random
import time
import datetime
#from scipy.stats import poisson


INFO_001 = "[INFO:001]"
INFO_002 = "[INFO:002]"
INFO_003 = "[INFO:003]"
INFO_004 = "[INFO:004]"
INFO_005 = "[INFO:005]"
INFO_006 = "[INFO:006]"
INFO_007 = "[INFO:007]"
INFO_008 = "[INFO:008]"
INFO_009 = "[INFO:009]"
INFO_010 = "[INFO:010]"
INFO_011 = "[INFO:011]"
INFO_012 = "[INFO:012]"


ERR_001 = "[ERR:001]"
ERR_002 = "[ERR:002]"
ERR_003 = "[ERR:003]"
ERR_004 = "[ERR:004]"
ERR_005 = "[ERR:005]"


def reverse_comp( seq ) :

    reversed = ''

    for base in "".join(seq)[::-1] : 
        if  base == 'A' : 
            reversed += 'T'  
        elif  base == 'C' : 
            reversed += 'G'  
        elif  base == 'T' : 
            reversed += 'A'  
        elif  base == 'G' : 
            reversed += 'C'  
        elif  base == 'N' : 
            reversed += 'A'  
        elif  base == 'a' : 
            reversed += 't'  
        elif  base == 'c' : 
            reversed += 'g'  
        elif  base == 'g' : 
            reversed += 'c'  
        elif  base == 't' : 
            reversed += 'a'  
        elif  base == 'n' : 
            reversed += 'a'  
        else : 
            print >> sys.stderr, "{0} not supported base : {1}".format( ERR_004, base )
   
    return reversed

    

def mutate( seq, seq_qual, sub_mu, in_mu, del_mu, rev_strd, replace ) :
    seq_list = list(seq)
    qual_list = list(seq_qual)

    i = 0
    #for i in range(len(seq)) :
    #print >> sys.stderr, "replace", replace

    if  sub_mu == 0 and in_mu == 0 and del_mu == 0 : 
        #print >> sys.stderr, "perfect"
        pass

    else :
        while True : 
    
            if  0 <= i < len(seq) :
                r = random.random()
        
                try : 
                    #print >> sys.stderr, i, r, sub_mu, in_mu, del_mu, seq[i]
         
                    if  r <= ( sub_mu ) : # substitution
                        #print random.randint(1,4) 
                        while True : 
                            base =  random.choice('ACGT') 
                            if  base != seq[i] : 
                                seq_list[i] = base
                                #qual_list[i] = 'B' 
                                qual_list[i] = '2' 
                                #print seq, "[sub]", len(seq)           
                                break
                            else : 
                                pass
                    elif  r <= ( sub_mu + in_mu ) : # insertion
                        base =  random.choice('ACGT') 
                        seq_list.insert(i, base )
                        seq = seq[:i] + base + seq[i:]
                        #qual_list.insert(i, '9' )
                        qual_list.insert(i, '1' )
                    elif  r <= ( sub_mu + in_mu + del_mu ) : # deletion
                        seq_list.pop(i)
                        seq = seq[:i] + seq[i+1:]
                        qual_list.pop(i)
                        #qual_list.insert(i,'7');
                        #print seq, "[del]", len(seq)
                    else : 
                        pass
        
                except IndexError : 
    	            print >> sys.stderr, "{0} i = {1} len(seq) = {2}".format( ERR_005, i, len(seq))
        
                i += 1
    
            else : 
                break
    
   
    if  rev_strd == 'on' :
        if  random.randint( 0, 1 ) == 0 : # forward
            if  replace == 'yes' : 
                return ( "F", "".join(seq_list).replace( 'N', 'A').replace( 'Y', 'A' ).replace( 'S', 'A'), "".join(qual_list) ) 
            else : 
                return ( "F", "".join(seq_list), "".join(qual_list) ) 
        else :  # backward, reversed strand
            reversed = reverse_comp("".join(seq_list))
            if  replace == 'yes' : 
                return ( "R", reversed.replace( 'N', 'A' ).replace( 'Y', 'A').replace( 'S', 'A' ), "".join(qual_list) ) 
            else :
                return ( "R", reversed, "".join(qual_list) ) 
    else : # forward 
        if  replace == 'yes' : 
            return ( "F", "".join(seq_list).replace( 'N', 'A' ).replace( 'Y', 'A').replace( 'S', 'A'), "".join(qual_list) ) 
        else :
            return ( "F", "".join(seq_list), "".join(qual_list) ) 


#def intro( fasta_file, read_mu, read_dist, cov_mu, cov_dist,
#           sub_mu, sub_dist, in_mu, in_dist, del_mu, del_dist ) :
def intro_readsim( tech, fasta_file, read_mu, read_dist, cov_mu,
           sub_mu, in_mu, del_mu ) :

    print >> sys.stderr, "================================================================================"
    print >> sys.stderr, "                         Read Simulator for Long Reads\n"
    print >> sys.stderr, "Sequencing Technology : {0}".format( tech )
    print >> sys.stderr, "FASTA file path : {0}".format( fasta_file )
    if  read_mu == 0 :
        pass
    else : 
        print >> sys.stderr, "Reads : mean({0}bp) distribution({1})".format( read_mu, read_dist )
    print >> sys.stderr, "Coverage : mean({0}x)".format( cov_mu )
    print >> sys.stderr, "Mutation (Substitution) : mean({0}%)".format( sub_mu*100 )
    print >> sys.stderr, "Mutation (Insertion)    : mean({0}%)".format(  in_mu*100 )
    print >> sys.stderr, "Mutation (Deletion)     : mean({0}%)".format( del_mu*100 )
    print >> sys.stderr, "================================================================================"
    

def output( outfile1, outfile2 ) : 

    print >> sys.stderr, "================================================================================"
    print >> sys.stderr, "{0} is generated.".format( outfile1 )
    print >> sys.stderr, "{0} is generated.".format( outfile2 )
    print >> sys.stderr, "================================================================================"


def get_num_from_exp( mu ) :

    return int( round( random.expovariate( 1.0/mu ) ) )

def get_num_from_normal( mu ) : 

    return int( round( random.gauss( mu, 1 ) ) )


def count_orientation( orientation, frd_cnt, bkd_cnt ) : 
    
    if  orientation == "F" : 
        frd_cnt += 1
    elif orientation == "R" : 
        bkd_cnt += 1
    else : 
        print >> sys.stderr, "[ERROR] Wrong orientation :", orientation 


    return frd_cnt, bkd_cnt

        




    
