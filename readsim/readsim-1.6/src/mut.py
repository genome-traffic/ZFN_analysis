#!/usr/bin/env python

# usage: ./cmd <fastafile> 

import sys
import random
import math
import common
from random import shuffle
import time
import datetime

MAX_LEN = 160000
MIN_LEN = 20


def mutate_write(  seq, read_len, read_num, pos, out_fa, out_fq, forward, backward,
                   sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, cov_mu,
                   sub_mu, in_mu, del_mu ) :

    if  (len(seq) - read_len) < 0 :
        read_len = len(seq)

    if  (len(seq) - read_len) >= 0 :
        if  sub_cmd == "fafq" :
            read_num += 1
            read_name = "{0}_{1:09d}_L{2:09d}:{3:09d}-{4:09d}".format(
                                                             prefix,
                                                             read_num, #str(read_num),
                                                             read_len,
                                                             pos,
                                                             pos + read_len - 1 )
            #print seq[0:read_len] 
            orientation, mut_seq, mut_qual =  common.mutate( seq[0:read_len], "M"*read_len, sub_mu, in_mu, del_mu, rev_strd, replace )
            forward, backward = common.count_orientation( orientation, forward, backward )

            out_fa.write( ">" + read_name + ":" + orientation + "\n" )
            out_fa.write( mut_seq + "\n" )
            out_fq.write( ">" + read_name + ":" + orientation + "\n" )
            out_fq.write( mut_qual + "\n" )

        elif  sub_cmd == "fq" :
            read_num += 1

            read_name = "{0}_{1:09d}_L{2:09d}:{3:09d}-{4:09d}".format(
                                                             prefix,
                                                             read_num, #str(read_num),
                                                             read_len,
                                                             pos,
                                                             pos + read_len - 1 )
            #print seq[0:read_len] 
            orientation, mut_seq, mut_qual =  common.mutate( seq[0:read_len], "M"*read_len, sub_mu, in_mu, del_mu, rev_strd, replace )
            forward, backward = common.count_orientation( orientation, forward, backward )
            out_fq.write( "@" + read_name + ":" + orientation + "\n" )
            out_fq.write( mut_seq + "\n" )
            out_fq.write( "+\n" )
            out_fq.write( mut_qual + "\n" )

        elif sub_cmd == "fa" :
            read_num += 1
            read_name = "{0}_{1:09d}_L{2:09d}:{3:09d}-{4:09d}".format(
                                                             prefix,
                                                             read_num, #str(read_num),
                                                             read_len,
                                                             pos,
                                                             pos + read_len - 1 )
            orientation, mut_seq, mut_qual =  common.mutate( seq[0:read_len], "M"*read_len, sub_mu, in_mu, del_mu, rev_strd, replace )
            forward, backward = common.count_orientation( orientation, forward, backward )
            out_fa.write( ">" + read_name + ":" + orientation + "\n" )
            out_fa.write( mut_seq + "\n" )
        else :
            pass

    else :
        #print >> sys.stderr, "{0} seq len = {1} < read len {2}".format( common.INFO_008, len(seq), read_len )
        pass

    return ( read_num, forward, backward )


def get_pos_len ( genlen, lens, c, times ) :

    read_len = lens[c]

    if  ( not read_len ) or ( read_len == 0 ) or ( read_len == "EOF" ) :
        c = 0
        read_len = lens[c]

    c += 1

    pos = random.randint( 0, genlen-1 )
    while pos + read_len > genlen :
        pos = random.randint( 0, genlen-1-read_len )

    read_len = int( read_len*times )  # here/change

    return ( pos, read_len, c )


def sim_reads_given_chr( sub_cmd, prefix, seq, lens, 
                         rev_strd, replace, read_mu, read_dist, cov_mu,
                         sub_mu, in_mu, del_mu, times ) :

    total_base = 0
    genlen = len(seq)
    i = 0
    c = 0
    read_num = 0 
    pos = 0 

    forward = 0
    backward = 0

    if  sub_cmd == "fq" :
        out_fq = open( prefix + ".fastq", 'w' )
        print >> sys.stderr, "{0} {1} {2} is created".format( common.INFO_002, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), prefix + ".fq" )
    elif  sub_cmd == "fa" :
        out_fa = open( prefix + ".fasta", 'w' )
        print >> sys.stderr, "{0} {1} {2} is created".format( common.INFO_003, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), prefix + ".fa" )
    elif  sub_cmd == "fafq" : # depreciated from v1.6
        out_fa = open( prefix + ".fa", 'w' )
        out_fq = open( prefix + ".fq", 'w' )
        print >> sys.stderr, "{0} {1} {2}, {3} are created".format( common.INFO_004, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), prefix + ".fa", prefix + ".fq" )
    else :
        print >> sys.stderr, "{0} sub command {1} is NOT supported".format( common.ERR_001, sub_cmd )


    if  read_dist == "uniform" : 
        read_len = read_mu 
        gap = int(read_mu/cov_mu)
        genlen = len(seq)
        init = 0

        while True :
            if  pos % gap == 1 :
                if  sub_cmd == "fq" :
                    read_num, forward, backward = mutate_write( seq[pos:pos+read_len], read_len, read_num, \
                                                                pos+init, out_fq, out_fq, forward, backward, \
                                                                sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, \
                                                                cov_mu, sub_mu, in_mu, del_mu )
                elif  sub_cmd == "fa" :
                    read_num, forward, backward = mutate_write( seq[pos:pos+read_len], read_len, read_num, \
                                                                pos+init, out_fa, out_fa, forward, backward, \
                                                                sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, \
                                                                cov_mu, sub_mu, in_mu, del_mu )
                elif  sub_cmd == "fafq" :
                    read_num, forward, backward = mutate_write( seq[pos:pos+read_len], read_len, read_num, \
                                                                pos+init, out_fa, out_fq, forward, backward, \
                                                                sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, \
                                                                cov_mu, sub_mu, in_mu, del_mu )
                else :
                    print >> sys.stderr, "{0} sub command {1} is NOT supported".format( common.ERR_001, sub_cmd )
                
                if  pos + read_len > genlen : 
                    break 
                else : 
                    #seq = seq[gap:]
                    pos += gap

            else :
                seq = seq[1:]
                pos += 1
                init = pos

            if  ( pos / 100000 ) >= i :
                i += 1
                print >> sys.stderr, "{0} {1} position {2} has been processed. ({3:.2f}x)".format( common.INFO_009, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), pos, float(pos)/genlen )
        print >> sys.stderr, "{0} {1} position {2} has been processed. ({3:.2f}x)".format( common.INFO_009, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), pos, float(pos)/genlen )

    else :  # normal/exp/read distribution files 

        while total_base < (cov_mu * genlen):

            pos, read_len, c = get_pos_len( genlen, lens, c, times ) 
            total_base += read_len

            #print >> sys.stderr, pos, read_len, total_base 

            if  sub_cmd == "fq" :
                read_num, forward, backward = mutate_write( seq[pos:pos+read_len], read_len, read_num, \
                                                            pos, out_fq, out_fq, forward, backward, \
                                                            sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, \
                                                            cov_mu, sub_mu, in_mu, del_mu )
            elif  sub_cmd == "fa" :
                read_num, forward, backward = mutate_write( seq[pos:pos+read_len], read_len, read_num,  \
                                                            pos, out_fa, out_fa, forward, backward, \
                                                            sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, \
                                                            cov_mu, sub_mu, in_mu, del_mu )
            elif  sub_cmd == "fafq" :
                read_num, forward, backward = mutate_write( seq[pos:pos+read_len], read_len, read_num, \
                                                            pos, out_fa, out_fq, forward, backward, \
                                                            sub_cmd, prefix, rev_strd, replace, read_mu, read_dist, \
                                                            cov_mu, sub_mu, in_mu, del_mu )
            else :
                print >> sys.stderr, "{0} sub command {1} is NOT supported".format( common.ERR_001, sub_cmd )

            #if  ( total_base / (genlen*100) ) >= i :
            if  ( total_base / 1000000 ) >= i :
                i += 1
                 #print >> sys.stderr, "{0} Total {1} bases has been generated.".format( common.INFO_009, total_base )
                print >> sys.stderr, "{0} {1} position {2} has been processed. ({3:.2f}x)".format( common.INFO_012, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), total_base, float(total_base)/genlen )
        print >> sys.stderr, "{0} {1} position {2} has been processed. ({3:.2f}x)".format( common.INFO_012, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), total_base, float(total_base)/genlen )


    if  sub_cmd == "fq" :
        out_fq.close()
    elif  sub_cmd == "fa" :
        out_fa.close()
    elif  sub_cmd == "fafq" :
        out_fa.close()
        out_fq.close()
    else :
        print >> sys.stderr, "{0} sub command {1} is NOT supported".format( common.ERR_002, sub_cmd )

    return ( forward, backward ) 

