#!/usr/bin/env python
import ast
import os
import sys
import re
import numpy as np
from tqdm import tqdm
import json


def fastq2fasta(f_dir,read_len,path):
    miseq_dir      = f_dir
    fragmented_dir = "%s/fragmented_fasta"%(path)

    fastq = get_fastq(miseq_dir)
    for F in fastq:
        make_fragmented_fasta("%s/%s"%(miseq_dir,F),500000,int(read_len),fragmented_dir) #File, split num
    pass

def exec_sh(sh_dir):
    fs = get_sh(sh_dir)
    for f in tqdm(fs):
        os.system('sh %s/%s'%(sh_dir,f))


def mk_blast_sh(fasta,db,seq_dir):
    n  =1
    q  =0
    fs = get_fasta(fasta)
    for f in fs:
        sh = seq_dir +'/blast/sh/BC_blast_%d.sh'%(n)
        with open(sh,"w") as F:
            name = f.split("/")[-1].split(".")[0]
            command = 'blastn -task blastn-short -strand plus -db  %s -outfmt 10 -evalue 1e-10 '%(db)
            blastsave = seq_dir+'/blast/out/'+name+'.blast'
            command += '-query '+'%s/%s'%(fasta,f)+' -out '+blastsave
            F.write(command)
        F.close()

        n+=1


def make_fragmented_fasta(F,split_num,read_length,out_dir):
    line_c  = 0
    read_c  = 0
    file_count   = 0
    with open(F,"r") as f:
        fasta_tup = []
        for line in f:
            #print(line)
            line = line.split("\n")[0]
            if (line_c % 4 == 0):
                read_ID = line.split(" ")[0].replace(":","_").replace("@","")

            elif (line_c % 4 == 1):
                if read_length > len(line):
                    seq = line
                else:
                    seq     = line[:read_length]

            elif (line_c % 4 == 3):
                qscore  = line
                fasta_tup.append( (read_ID,seq) )
                read_c +=1

                if (read_c % split_num == 0):
                    file_count +=1


                    f_name = ("").join([F.split("/")[-1].split(".")[0] , "_%d"%(file_count),".fna"])
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
            out.append(("c%s-DNTAG"%(BC_ID),c_dntag))

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
            out.append(("c%s-DNTAG"%(BC_ID),c_dntag))
    make_fasta_from_listoftups(out,"%s/%s_database.fna"%(db_dir,run_name))

    return "%s/%s_database.fna"%(db_dir,run_name)

def csv2LL(f_name):
    LL = []
    with open(f_name,"r") as F:

        for line in F:
            cols = line.split("\n")[0].split(",")
            LL.append(cols)
    return LL

def get_fastq(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-6:] == ".fastq" :
            files.append(i)
    return files

def get_fasta(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-4:] == ".fna" :
            files.append(i)
    return files


def get_sh(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-3:] == ".sh" :
            files.append(i)
    return files



def get_blast(x):
    files = []
    lis = os.listdir("%s" %x)
    for i in lis:
        if i[-6:] == ".blast" :
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
            f.write('>%s\n'%(str(I[0])))
            f.write('%s\n'%(str(I[1])))
        f.close()
        print("Made fna file : %s"%(name))




def parse_blast_out(PATH,run_name,db_fna):
    seq_dir = "%s/Data/%s"%(PATH,run_name)
    BC_fasta      = db_fna
    bar2num_file  = '%s/bar2num.csv'%(PATH)


    seq_DB  = fasta2dict(BC_fasta)
    #print(seq_DB)
    bar2num = bar2num_d(bar2num_file)

    #print Dumper $seq_DB;
    list_of_fs   = get_blast('%s/blast/out'%(seq_dir))
    list_of_f_id = [f.split(".blast")[0].replace("R2","RX") for f in list_of_fs if ("R2" in f)]
    D = {}
    seq_reads = {}
    seq_reads["R1"] = {}
    seq_reads["R2"] = {}


    for f in tqdm(list_of_f_id):
        BC_R1_blast = seq_dir+'/blast/out/'       +f.replace("RX","R1")+'.blast'
        BC_R2_blast = seq_dir+'/blast/out/'       +f.replace("RX","R2")+'.blast'
        fasta_R1    = seq_dir+'/fragmented_fasta/'+f.replace("RX","R1")+'.fna'
        fasta_R2    = seq_dir+'/fragmented_fasta/'+f.replace("RX","R2")+'.fna'

        seq_reads["R1"].update(fasta2dict(fasta_R1))
        seq_reads["R2"].update(fasta2dict(fasta_R2))


        files = [BC_R1_blast,BC_R2_blast]


        for f in files:
            if "R1" in f:
                read_direction = "R1"
            elif "R2" in f:
                read_direction = "R2"


            with open(f,"r") as F:
                for line in F:

                    cols        = line.split(",")
                    read        =  cols[0]
                    template    =    cols[1]
                    bit         =    float(cols[2])
                    str_on_read =    int(cols[6])
                    end_on_read =    int(cols[7])
                    str         =    int(cols[8])
                    end         =    int(cols[9])
                    evalue      =    float(cols[10])
                    if( (str_on_read<end_on_read )and (bit>=90) and (str <= 20 ) and ( end >= len(seq_DB[template])-20)):
                        #print(read_direction)
                        el = {}
                        el["str"]      = str_on_read - str +1
                        el["end"]      = end_on_read + len(seq_DB[template]) - end
                        el["evalue"]   = evalue
                        el["bitscore"] = bit
                        if("TAG" in template):
                            el["name"] = template

                        try:
                            D[read][read_direction][template] = el
                        except KeyError:
                            try:
                                D[read][read_direction] = {}
                                D[read][read_direction][template] = el
                            except KeyError:
                                D[read] = {}
                                D[read][read_direction] = {}
                                D[read][read_direction][template] = el
                        #print(template,D[read][read_direction][template])

            F.close()

    m = []
    for read in D:
        if len(m)<10:
            m.append(read)
        for read_direction in D[read]:
            data = D[read][read_direction]
            ranked = [(template,data[template]["bitscore"]) for template in data]
            ranked.sort(key=lambda x:x[1],reverse=True)
            D[read][read_direction] = data[ranked[0][0]]

    #print([D[i] for i in m])
    D2 = {}
    for read in D:
        if "R1" in D[read]:
            if "R2" in D[read]:
                R1       = D[read]["R1"]
                R2       = D[read]["R2"]

                seq_R1   = seq_reads["R1"][read]
                seq_R2   = seq_reads["R2"][read]

                if R1["str"] <=20:
                    if R2["str"] <=20:
                        D2[read] = {}
                        D2[read]["R1"] = {}
                        D2[read]["R2"] = {}
                        name   = R1["name"]
                        D2[read]["R1"]["name"]   = name
                        P1_TAG = seq_R1[:R1["str"]-1]
                        D2[read]["R1"]["P1_TAG"] = P1_TAG
                        name   = R2["name"]
                        D2[read]["R2"]["name"]   = name
                        P2_TAG = seq_R2[:R2["str"]-1]
                        D2[read]["R2"]["P2_TAG"] = P2_TAG

    count = {}
    print("Merging barcode count data from all BLAST outputs....")
    for read in tqdm(D2):
        P1_num = 0
        P2_num = 0
        P1_TAG = D2[read]["R1"]["P1_TAG"][-9:]
        P1_num = int(barcode_matching(bar2num,P1_TAG))
        P2_TAG = D2[read]["R2"]["P2_TAG"][-9:]
        P2_num = int(barcode_matching(bar2num,P2_TAG))


        if P1_num !=0:
            if P2_num !=0:
                if "name" in D2[read]["R1"]:
                    if "name" in D2[read]["R2"]:
                        P_TAG_comb = "P%02d-P%02d"%(P1_num,P2_num)
                        if D2[read]["R1"]["name"][0] =="c":
                            name1    = D2[read]["R1"]["name"][1:]
                        else:
                            name1    = D2[read]["R1"]["name"]
                        name1_BC = name1.split("-")[0]
                        name1_bfg= name1.split("-")[-1]
                        if D2[read]["R2"]["name"][0] =="c":
                            name2    = D2[read]["R2"]["name"][1:]
                        else:
                            name2    = D2[read]["R2"]["name"]
                        name2_BC = name2.split("-")[0]
                        name2_bfg= name2.split("-")[-1]

                        if(name1_bfg == name2_bfg):
                            fuse_type = 'NaN'
                            if name1_bfg == 'UPTAG':
                                fuse_type = 'UpUp'
                            if name2_bfg == 'DNTAG':
                                fuse_type = 'DnDn'

                        try:
                            count[P_TAG_comb][name1_BC][name2_BC][fuse_type] +=1
                        except KeyError:
                            try:
                                count[P_TAG_comb][name1_BC][name2_BC][fuse_type] = 1
                            except KeyError:
                                try:
                                    count[P_TAG_comb][name1_BC][name2_BC] = {}
                                    count[P_TAG_comb][name1_BC][name2_BC][fuse_type] = 1
                                except KeyError:
                                    try:
                                        count[P_TAG_comb][name1_BC] = {}
                                        count[P_TAG_comb][name1_BC][name2_BC] = {}
                                        count[P_TAG_comb][name1_BC][name2_BC][fuse_type] = 1
                                    except KeyError:
                                            count[P_TAG_comb] = {}
                                            count[P_TAG_comb][name1_BC] = {}
                                            count[P_TAG_comb][name1_BC][name2_BC] = {}
                                            count[P_TAG_comb][name1_BC][name2_BC][fuse_type] = 1
    print("Outputing the count data to %s/Data/%s/barcode_counts/counts.txt"%(PATH,run_name))
    with open('%s/Data/%s/barcode_counts/counts.txt'%(PATH,run_name), 'w') as file:
        file.write(json.dumps(count))
    print("\nDone")


def fasta2dict(f):
    d = {}
    with  open(f,"r") as F:
        for line in F:
            if line[0]==">":
                name = line.split(">")[1].split("\n")[0]
            else:
                #print(name)
                d[name] = line.split("\n")[0]
    F.close()
    return d




def barcode_matching(bar2num,seq):
    bar_d = {}
    count = {}
    for k  in bar2num:
        count = {}
        for i in range(len(k)):
            k_nuc = k[i]
            for j in range(len(seq)):
                seq_nuc = seq[j]
                rel = i-j
                if(k_nuc == seq_nuc):
                    try:
                        count[rel] +=1
                    except KeyError:
                        count[rel] =1

        if len(count) >0:
            largest_hit = max([count[rel] for rel in count])
            bar_d[k] = largest_hit
        else:
            bar_d[k] = 0

    val = 0

    bar_counts = [(k,bar_d[k]) for k in bar_d if (bar_d[k]>5)]
    bar_counts.sort(key=lambda x:x[1],reverse=True)
    #print(bar_counts)
    if len(bar_counts) >1:
        if bar_counts[0][1] == bar_counts[1][1]:
            val = 0
        else:
            val = bar2num[bar_counts[0][0]]
    elif len(bar_counts) ==1:
        val = bar2num[bar_counts[0][0]]
    else:
        val = 0

    return val

def bar2num_d(f):
    d = {}
    with open(f,"r") as F:
        for line in F:
            cols = line.split("\n")[0].split(",")
            d[cols[1]] = cols[0]
    F.close()
    return d
