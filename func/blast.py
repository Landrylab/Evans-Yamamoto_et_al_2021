#!/usr/bin/env python
import ast
import os
import sys
import re
import numpy as np
from tqdm import tqdm


def fastq2fasta(f_dir,read_len,path):
    miseq_dir      = f_dir
    fragmented_dir = "%s/fragmented_fasta"%(path)

    fastq = get_files(miseq_dir)
    for F in fastq:
        make_fragmented_fasta("%s/%s"%(miseq_dir,F),100000,int(read_length),fragmented_dir) #File, split num
    pass


def make_fragmented_fasta(F,split_num,read_length,out_dir):
    line_c  = 0
    read_c  = 0
    file_count   = 0
    with open(F,"r") as f:
        fasta_tup = []
        for line in f:
            #print(line)
            if (line_c % 4 == 0):
                read_ID = line.replace(" ","_").replace(":","_").replace("@","")

            elif (line_c % 4 == 1):
                if read_length > line:
                    seq = line
                else:
                    seq     = line[:read_length]

            elif (line_c % 4 == 3):
                qscore  = line
                fasta_tup.append( (read_ID,seq) )
                read_c +=1

                if (read_c % split_num == 0):
                    file_count +=1


                    f_name = ("").join([F.split("/")[-1].split(".")[0] , "_6%d"%(file_count),".fna"])
                    out    = ("/").join( [out_dir,f_name])

                    make_fasta_from_listoftups(fasta_tup, out )

                    fasta_tup = []
            line_c += 1
        f.close()
    pass


def mkdb(db_file,const_seq,run_name,db_dir):
    const_data = csv2LL(const_seq)
    #print(const_data)
    const_d = {}
    for i in const_data:
        const_d[i[0]] = i[1]
        const_d["c%s"%(i[0])] = rev_comp(i[1])

    DB =  csv2LL(db_file)

    out = []

    for i in DB:
        type  = i[0]
        BC_ID = i[5]
        BC1   = i[6]
        BC2   = i[7]

        if(type == "DHFR3" ):
            type = "DB"

            U1 = type + 'U1-primer'
            U2 = 'c' + type + 'U2-primer'
            D1 = type + 'D1-primer'
            D2 = 'c' + type + 'D2-primer'

            uptag = const_d[U1] + BC1 + const_d[U2]
            dntag = const_d[D1] + BC2 + const_d[D2]

            c_uptag = rev_comp(uptag)

            c_dntag = rev_comp(dntag)

            out.append(("%s-UPTAG"%(BC_ID),uptag))
            out.append(("c%s-UPTAG"%(BC_ID),c_uptag))

            out.append(("%s-DNTAG"%(BC_ID),dntag))
            out.append(("%s-DNTAG"%(BC_ID),c_dntag))

        if(type == "DHFR12" ):
            type = "AD"

            U1 = type + 'U1-primer'
            U2 = 'c' + type + 'U2-primer'
            D1 = type + 'D1-primer'
            D2 = 'c' + type + 'D2-primer'


            uptag = const_d[U1] + BC1 + const_d[U2]
            dntag = const_d[D1] + BC2 + const_d[D2]

            c_uptag = rev_comp(uptag)
            c_dntag = rev_comp(dntag)

            out.append(("%s-UPTAG"%(BC_ID),uptag))
            out.append(("c%s-UPTAG"%(BC_ID),c_uptag))

            out.append(("%s-DNTAG"%(BC_ID),dntag))
            out.append(("%s-DNTAG"%(BC_ID),c_dntag))
    make_fasta_from_listoftups(out,"%s/%s_database.fna"%(db_dir,run_name))

    return "%s/%s_database.fna"%(db_dir,run_name)

def csv2LL(f_name):
    LL = []
    with open(f_name,"r") as F:

        for line in F:
            cols = line.split("\n")[0].split(",")
            LL.append(cols)
    return LL

def get_files(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-6:] == ".fastq" :
            files.append(i)
    return files


def rev_comp(seq):
    complement_dict = {'A':'T','T':'A','G':'C','C':'G',"N":"N"}
    return "".join([complement_dict[base] for base in reversed(seq)])


def generate_directory(mypath):
    if not os.path.isdir(mypath):
        os.makedirs(mypath)


def make_fasta_from_listoftups(L,name):
    with open("%s"%name,"w") as f:
        for I in L:
            #print I
            f.write('>%s'%(str(I[0])))
            f.write('%s\n'%(str(I[1])))
        f.close()
        print("Made fna file : %s"%(name))
