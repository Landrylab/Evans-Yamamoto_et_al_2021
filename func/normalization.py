# coding: utf-8
import ast
import os
import sys
import re
from pprint import pprint
import numpy as np
import itertools


#############################
## Defining small functions
#############################

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r
def csv2LL(f_name):
    LL = [[]]
    with open(f_name,"r") as F:
        c =1
        for line in F:
            cols = line.split("\n")[0].split(",")
            if (c> 1):
                LL.append(cols)
            c+=1
    return LL
def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    # This script outputs dict as string, so do not print anything
    print("Generated : %s"%(name))
    F.close()
def reading(name):
    with open(name, 'r') as F:
        ob_from_file = F.read().replace('\n','')
        mydata = ast.literal_eval(ob_from_file)
    F.close()
    return mydata

def get_sample_name(nm,tag_dic):
    for tag in tag_dic:
        if (tag_dic[tag]["Control"] == "%s*"%(nm)):
            sample = tag_dic[tag]["Sample"]
    return sample

def get_mapdata(f):
    D = {}

    with open(f,"r") as F:
        for line in F:
            cols = line.split(",")
            d = {
                    "Type"        : cols[0],
                    "Category"    : cols[1],
                    "ID"          : cols[2],
                    "ORF_SYMBOL"  : cols[3],
                    "BC_ver"      : cols[4],
                    "BC_ID"       : cols[5]
                }
            try:
                D[cols[0]].append(d)
            except KeyError:
                D[cols[0]] = [d]
        F.close()
    return D

def get_tagdata(f):
    D = {}
    with open(f,"r") as F:
        for line in F:
            cols = line.split(",")
            D[cols[0]] = {
                            "Screen"      : cols[1],
                            "Sample"      : cols[2],
                            "Timepoint"   : cols[3],
                            "Comment"     : cols[4],
                            "Control"     : cols[5].split("\n")[0]
                        }
        F.close()
    return D

def tot_reads(counts,multiplex_tag):
    indecies = [i for i in multiplex_tag]
    indecies.sort()
    for tag in indecies:
        c = 0
        st_n = 0
        for st in counts[count][tag]:
            for st2 in counts[count][tag][st]:
                for bfg in counts[count][tag][st][st2]:
                    c += counts[count][tag][st][st2][bfg]
        print(tag,multiplex_tag[tag]["Sample"],"Total reads %d "%(c))

def output_raw_reads(counts,tag,norm_dir):
    k = list(tag.keys())
    k.sort()
    header = []
    for iu in k:
        for bfg in ["UpUp","DnDn"]:
            header.append((".").join([iu,bfg]))

    LL = [["Bait_BC_ID","Prey_BC_ID"]]
    LL[0]+= header
    D = {}
    for count_dat in k:
        for bait in counts[count_dat]:
            for prey in counts[count_dat][bait]:
                for bfg in counts[count_dat][bait][prey]:
                    try:
                        D[(bait,prey)][count_dat][bfg] = counts[method][count_dat][bait][prey][bfg]
                    except KeyError:
                        try:
                            D[(bait,prey)][count_dat] = {}
                            D[(bait,prey)][count_dat][bfg] = counts[method][count_dat][bait][prey][bfg]
                        except KeyError:
                            D[(bait,prey)] = {}
                            D[(bait,prey)][count_dat] = {}
                            D[(bait,prey)][count_dat][bfg] = counts[method][count_dat][bait][prey][bfg]




    for i in D:
        l = [i[0],i[1]]
        for ind in k:
            for bfg in ["UpUp","DnDn"]:
                #print(ind,bfg)
                try:
                    val = D[i][ind][bfg]
                except KeyError:
                    val = 0
                l.append(val)
        LL.append(l)
    LL2csv(LL,"%s/raw_counts.csv"%(norm_dir))

#########################
## Function to compute organize dict and sum strain abundance
#########################

def organize_data(data,map_table,tag_data,method):
    data2 = {}
    sums = {}
    alpha  = 1
    barcode_fusion_type = ["UpUp","DnDn"]
    for index in tag_data:
        if index in data.keys():
            selection = ("_").join([tag_data[index]["Sample"],tag_data[index]["Screen"]])
            control   = tag_data[index]["Control"]
            if control[-1] == "*":
                selection = control[:-1]
            print(index,selection)
            data2[selection] = {}

            if method== "Y2H":
                for bait in map_table["DB"]:
                    bait_ID = bait["ID"]
                    bait_BC_ID = bait["BC_ID"]
                    bait_BC_ver = bait["BC_ver"]
                    if bait_ID not in data2[selection].keys():
                        data2[selection][bait_ID] = {}

                    for prey in map_table["AD"]:
                        prey_ID = prey["ID"]
                        prey_BC_ID = prey["BC_ID"]
                        prey_BC_ver = prey["BC_ver"]
                        if prey_ID not in data2[selection][bait_ID].keys():
                            data2[selection][bait_ID][prey_ID] = {}

                        for bfg in barcode_fusion_type:

                            try:
                                raw   = data[index][prey_BC_ID][bait_BC_ID][bfg]
                            except KeyError:
                                raw = 0

                            count = raw + alpha

                            try:
                                data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg] = {}
                                data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["raw"] = raw
                                data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"] = count
                            except KeyError:
                                try:
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["raw"] = raw
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"] = count

                                except KeyError:
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["raw"] = raw
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"] = count

            if method== "PCA":
                for bait in map_table["DHFR3"]:
                    bait_ID = bait["ID"]
                    bait_BC_ID = bait["BC_ID"]
                    bait_BC_ver = bait["BC_ver"]
                    if bait_ID not in data2[selection].keys():
                        data2[selection][bait_ID] = {}

                    for prey in map_table["DHFR12"]:
                        prey_ID = prey["ID"]
                        prey_BC_ID = prey["BC_ID"]
                        prey_BC_ver = prey["BC_ver"]
                        if prey_ID not in data2[selection][bait_ID].keys():
                            data2[selection][bait_ID][prey_ID] = {}

                        for bfg in ["UpUp","DnDn"]:

                            try:
                                raw   = data[index][prey_BC_ID][bait_BC_ID][bfg]

                            except KeyError:
                                raw = 0
                            count = raw + alpha

                            try:
                                data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg] = {}
                                data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["raw"] = raw
                                data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"] = count
                            except KeyError:
                                try:
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["raw"] = raw
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"] = count

                                except KeyError:
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg] = {}
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["raw"] = raw
                                    data2[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"] = count

    return data2

#############################
## Function to remove strains absent in >90% of diploids
#############################

def absent(data,absent_th,method):
    absent_rate = {}
    raw_count   = {}

    for selection in data:
        if selection.split("_")[0] == "Control":
            for bait in data[selection]:
                for prey in data[selection][bait]:
                    for bait_bc_ver in data[selection][bait][prey]:
                        raw_count["%s_%s"%(bait,bait_bc_ver)] = {}
                        for prey_bc_ver in data[selection][bait][prey][bait_bc_ver]:
                            for bfg in ["UpUp","DnDn"]:
                                try:
                                    raw_count["%s_%s"%(bait,bait_bc_ver)][bfg].append( data[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["raw"] )
                                except KeyError:
                                    raw_count["%s_%s"%(bait,bait_bc_ver)][bfg] = []
                                    raw_count["%s_%s"%(bait,bait_bc_ver)][bfg].append(data[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["raw"])
                                try:
                                    raw_count["%s_%s"%(prey,prey_bc_ver)][bfg].append(data[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["raw"])
                                except KeyError:
                                    try:
                                        raw_count["%s_%s"%(prey,prey_bc_ver)][bfg] = []
                                        raw_count["%s_%s"%(prey,prey_bc_ver)][bfg].append(data[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["raw"])
                                    except KeyError:
                                        raw_count["%s_%s"%(prey,prey_bc_ver)] = {}
                                        raw_count["%s_%s"%(prey,prey_bc_ver)][bfg] = []
                                        raw_count["%s_%s"%(prey,prey_bc_ver)][bfg].append(data[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["raw"])

                                try:
                                    absent_rate["%s_%s"%(bait,bait_bc_ver)][bfg] = float(raw_count["%s_%s"%(bait,bait_bc_ver)][bfg].count(0))/len(raw_count["%s_%s"%(bait,bait_bc_ver)][bfg])
                                except KeyError:
                                    absent_rate["%s_%s"%(bait,bait_bc_ver)] = {}
                                    absent_rate["%s_%s"%(bait,bait_bc_ver)][bfg] = float(raw_count["%s_%s"%(bait,bait_bc_ver)][bfg].count(0))/len(raw_count["%s_%s"%(bait,bait_bc_ver)][bfg])
                                try:
                                    absent_rate["%s_%s"%(prey,prey_bc_ver)][bfg] = float(raw_count["%s_%s"%(prey,prey_bc_ver)][bfg].count(0))/len(raw_count["%s_%s"%(prey,prey_bc_ver)][bfg])
                                except KeyError:
                                    absent_rate["%s_%s"%(prey,prey_bc_ver)] = {}
                                    absent_rate["%s_%s"%(prey,prey_bc_ver)][bfg] = float(raw_count["%s_%s"%(prey,prey_bc_ver)][bfg].count(0))/len(raw_count["%s_%s"%(prey,prey_bc_ver)][bfg])
            absent = {
            "bait":[],#(bait,bc_ver)
            "prey":[]
            }
            #pprint(selection,absent_rate)
            for haploid in absent_rate:
                ave  = (absent_rate[haploid]["UpUp"] + absent_rate[haploid]["DnDn"]) /2
                #print haploid,ave
                if method == "PCA":
                    if ave > absent_th:
                        if "DHFR12" in haploid:
                            #print haploid,"DHFR12"
                            absent["prey"].append(haploid)
                        if "DHFR3" in haploid:
                            absent["bait"].append(haploid)
                            #print haploid,"DHFR3"
                if method == "Y2H":
                    if ave > absent_th:
                        if "AD" in haploid:
                            #print haploid,"DHFR12"
                            absent["prey"].append(haploid)
                        if "DB" in haploid:
                            absent["bait"].append(haploid)
                            #print haploid,"DHFR3"




            for bait in data[selection]:
                for prey in data[selection][bait]:
                    for bait_bc_ver in data[selection][bait][prey]:
                        haploid = "%s_%s"%(bait,bait_bc_ver)
                        if haploid in absent["bait"]:
                            data[selection][bait][prey] = removekey(data[selection][bait][prey],bait_bc_ver)

                        else:
                            for prey_bc_ver in data[selection][bait][prey][bait_bc_ver]:
                                haploid = "%s_%s"%(prey,prey_bc_ver)
                                if haploid in absent["prey"]:
                                    data[selection][bait][prey][bait_bc_ver] = removekey(data[selection][bait][prey][bait_bc_ver],prey_bc_ver)

                            if len(data[selection][bait][prey][bait_bc_ver]) == 0:
                                data[selection][bait][prey] = removekey(data[selection][bait][prey],bait_bc_ver)

                    if len(data[selection][bait][prey]) == 0:
                        data[selection][bait] = removekey(data[selection][bait],prey)

                if len(data[selection][bait]) == 0:
                    data[selection] = removekey(data[selection],bait)


    return data


def get_sums(data,tag_data,method,map_table):
    sums = {}

    for selection in data:

        sums[selection] = {}
        sums[selection]["UpUp"] = {}
        sums[selection]["DnDn"] = {}
        sums[selection]["UpUp"]["all"] = 0
        sums[selection]["DnDn"]["all"] = 0

        if method== "Y2H":
            for bait in map_table["DB"]:
                bait_ID = bait["ID"]
                bait_BC_ID = bait["BC_ID"]
                bait_BC_ver = bait["BC_ver"]

                for prey in map_table["AD"]:
                    prey_ID = prey["ID"]
                    prey_BC_ID = prey["BC_ID"]
                    prey_BC_ver = prey["BC_ver"]

                    for bfg in ["UpUp","DnDn"]:

                        try:
                            count = data[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"]
                            try:
                                sums[selection][bfg][bait_ID][bait_BC_ver]  += count
                            except KeyError:
                                try:
                                    sums[selection][bfg][bait_ID][bait_BC_ver]  = count
                                except KeyError:
                                    sums[selection][bfg][bait_ID]  = {}
                                    sums[selection][bfg][bait_ID][bait_BC_ver]  = count

                            try:
                                sums[selection][bfg][prey_ID][prey_BC_ver]  += count
                            except KeyError:
                                try:
                                    sums[selection][bfg][prey_ID][prey_BC_ver]  = count
                                except KeyError:
                                    sums[selection][bfg][prey_ID]  = {}
                                    sums[selection][bfg][prey_ID][prey_BC_ver]  = count
                            sums[selection][bfg]["all"] += count
                        except KeyError:
                            pass


        if method== "PCA":
            for bait in map_table["DHFR3"]:
                bait_ID = bait["ID"]
                bait_BC_ID = bait["BC_ID"]
                bait_BC_ver = bait["BC_ver"]

                for prey in map_table["DHFR12"]:
                    prey_ID = prey["ID"]
                    prey_BC_ID = prey["BC_ID"]
                    prey_BC_ver = prey["BC_ver"]

                    for bfg in ["UpUp","DnDn"]:
                        try:
                            count = data[selection][bait_ID][prey_ID][bait_BC_ver][prey_BC_ver][bfg]["count"]
                            try:
                                sums[selection][bfg][bait_ID][bait_BC_ver]  += count
                            except KeyError:
                                try:
                                    sums[selection][bfg][bait_ID][bait_BC_ver]  = count
                                except KeyError:
                                    sums[selection][bfg][bait_ID]  = {}
                                    sums[selection][bfg][bait_ID][bait_BC_ver]  = count

                            try:
                                sums[selection][bfg][prey_ID][prey_BC_ver]  += count
                            except KeyError:
                                try:
                                    sums[selection][bfg][prey_ID][prey_BC_ver]  = count
                                except KeyError:
                                    sums[selection][bfg][prey_ID]  = {}
                                    sums[selection][bfg][prey_ID][prey_BC_ver]  = count

                            sums[selection][bfg]["all"] += count
                        except KeyError:
                            pass




    return sums

##############################
## Function to normalize scores based on strain abundance of bait and prey haploids in Control condition
##############################

def compute_s(data2,sums,map_table,tag_data,method):
    hap = {}
    data3 = {}

    for selection in data2:
        data3[selection] = {}
        if selection.split("_")[0] != "Control":
            control  = [ tag_data[index]["Control"] for index in tag_data if (("_").join([tag_data[index]["Sample"],tag_data[index]["Screen"]])==selection) ][0]

            for bait in  data2[control]:
                data3[selection][bait] = {}
                for prey in data2[control][bait]:
                    data3[selection][bait][prey] = {}
                    for bait_bc_ver in data2[control][bait][prey]:
                        data3[selection][bait][prey][bait_bc_ver] = {}
                        for prey_bc_ver in data2[control][bait][prey][bait_bc_ver]:
                            data3[selection][bait][prey][bait_bc_ver][prey_bc_ver] = {}


                            for bfg in ["UpUp","DnDn"]:

                                count = float(data2[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["count"])
                                total_count = float(sums[selection][bfg]["all"])
                                #i is bait, j is prey

                                F_ij = count/total_count


                                count_i = float(sums[control][bfg][bait][bait_bc_ver])
                                F_est_i     = count_i/total_count

                                count_j = float(sums[control][bfg][prey][prey_bc_ver])
                                F_est_j     = count_j/total_count

                                F_est_ij = F_est_i * F_est_j
                                s        = F_ij / F_est_ij






                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg] = data2[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]
                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["F_est"] = F_est_ij
                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["F"] = F_ij
                                if method =="PCA":
                                    data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["s"] = s+1
                                if method =="Y2H":
                                    data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["s"] = s
                                #pprint(data3[selection][bait][prey][bait_bc_ver][prey_bc_ver])

        else:
            hap[selection] = {}
            hap[selection]["bait"] = {}
            hap[selection]["prey"] = {}

            for bait in  data2[selection]:
                data3[selection][bait] = {}
                for prey in data2[selection][bait]:
                    data3[selection][bait][prey] = {}
                    for bait_bc_ver in data2[selection][bait][prey]:
                        data3[selection][bait][prey][bait_bc_ver] = {}
                        for prey_bc_ver in data2[selection][bait][prey][bait_bc_ver]:
                            data3[selection][bait][prey][bait_bc_ver][prey_bc_ver] = {}
                            FI = []
                            FJ = []
                            for bfg in ["UpUp","DnDn"]:
                                count = float(data2[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["count"])
                                total_count = float(sums[selection][bfg]["all"])
                                #i is bait, j is prey

                                count = float(data2[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["count"])
                                total_count = float(sums[selection][bfg]["all"])

                                F_ij = count/total_count


                                count_i = float(sums[selection][bfg][bait][bait_bc_ver])
                                F_est_i     = count_i/total_count

                                count_j = float(sums[selection][bfg][prey][prey_bc_ver])
                                F_est_j     = count_j/total_count

                                F_est_ij = F_est_i * F_est_j
                                hap[selection]["bait"]["%s-%s"%(bait,bait_bc_ver)] = F_est_i
                                hap[selection]["prey"]["%s-%s"%(prey,prey_bc_ver)] = F_est_j


                                FI.append(F_est_i)
                                FJ.append(F_est_j)
                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg] = data2[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]
                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["F_est"] = F_est_ij
                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["F"] = F_ij

                            Fi = np.average(FI)
                            Fj = np.average(FJ)

                            hap[selection]["bait"]["%s-%s"%(bait,bait_bc_ver)] = Fi
                            hap[selection]["prey"]["%s-%s"%(prey,prey_bc_ver)] = Fj



    return data3,hap

def count_haploids(PCA_hap,out_dir):

    HAP = [["Screeening","Ori","Strain","F"]]

    reps = {}

    for sel in PCA_hap:
        reps[sel] = {"bait":{},"prey":{}}

        b = 0
        p = 0
        b_orf = []
        p_orf = []
        for strain in PCA_hap[sel]["bait"]:
            l = ["PCA %s"%(sel),"Bait",strain,PCA_hap[sel]["bait"][strain]]
            HAP.append(l)
            b += 1
            #print (strain)
            bait = list(set([i["ORF_SYMBOL"] for i in PCA_db[strain.split("_")[0]] if (i['ID']==strain.split("-")[0])] ))[0]
            b_orf.append(bait)

            try:
                reps[sel]["bait"][bait].append(1)
            except KeyError:
                reps[sel]["bait"][bait] = [1]

        for strain in PCA_hap[sel]["prey"]:
            l = ["PCA %s"%(sel),"Prey",strain,PCA_hap[sel]["prey"][strain]]
            HAP.append(l)
            p +=1
            prey = list(set([i["ORF_SYMBOL"] for i in PCA_db[strain.split("_")[0]] if (i['ID']==strain.split("-")[0])] ))[0]
            p_orf.append(prey)
            try:
                reps[sel]["prey"][prey].append(1)
            except KeyError:
                reps[sel]["prey"][prey] = [1]

        print(sel,"\t\tBait: , (%s ORFs)"%(len(set(b_orf))),b,"Prey, (%s ORFs)"%(len(set(p_orf))),p)
    LL2csv(HAP,"%s/Haploid_F.csv"%(out_dir))
    return reps

def haploid_replicates(reps,out_dir):
    counts= [["Method","Cond","Rep","ori","n"]]

    for sel in reps:
        for ori in reps[sel]:
            for strain in reps[sel][ori]:
                rep = sel.split("_")[1]
                con = sel.split("_")[0]
                l = [con,rep,ori,len(reps[sel][ori][strain])]
                counts.append(l)
    LL2csv(counts,"%s/haplliod_n.csv"%(out_dir))

################################
## Normalize signal 's'
## For Y2H ; Normalize data by Nth quantile of bait 'median distribution' of 's' score
## For PCA ; Normalize data by median 's' score of both bait and prey


def compute_ds(data3,sums,map_table,tag_data,method):
    #Retrieving scores from matrix
    haploid_s_scores = {}
    for selection in data3:
        if selection.split("_")[0] != "Control":
            haploid_s_scores[selection] = {}
            #Estimating strain abundance based on control condition
            for bait in  data3[selection]:
                for prey in data3[selection][bait]:
                    for bait_bc_ver in data3[selection][bait][prey]:
                        for prey_bc_ver in data3[selection][bait][prey][bait_bc_ver]:
                            for bfg in data3[selection][bait][prey][bait_bc_ver][prey_bc_ver]:
                                s_ij = data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["s"]
                                #print selection,bait,prey,bait_bc_ver,prey_bc_ver,bfg,s_ij

                                try:
                                    haploid_s_scores[selection]["Prey_%s_%s"%(prey,prey_bc_ver)]["scores"].append(s_ij)
                                except KeyError:
                                    haploid_s_scores[selection]["Prey_%s_%s"%(prey,prey_bc_ver)] = {
                                    "scores": [s_ij]
                                    }

                                try:
                                    haploid_s_scores[selection]["Bait_%s_%s"%(bait,bait_bc_ver)]["scores"].append(s_ij)
                                except KeyError:
                                    haploid_s_scores[selection]["Bait_%s_%s"%(bait,bait_bc_ver)] = {
                                    "scores": [s_ij]
                                    }

    #Calculating median of scores for erach haploid strain, and computing beta based on supplied beta threashold.
    betas = {}
    for selection in haploid_s_scores:

        #print selection
        betas[selection]= {"Bait":{},"Prey":{}}
        for haploid in haploid_s_scores[selection]:
            strain_type = haploid.split("_")[0]
            betas[selection][strain_type][haploid] = {}
            haploid_s_scores[selection][haploid]["beta"] = {}
            haploid_s_scores[selection][haploid]["median"] = np.median( haploid_s_scores[selection][haploid]["scores"] )
            #print haploid, np.median( haploid_s_scores[selection][haploid]["scores"] )
            haploid_s_scores[selection][haploid]["s_minus_median"] = sorted(list(haploid_s_scores[selection][haploid]["scores"] - haploid_s_scores[selection][haploid]["median"]),reverse =True)
            positive_values = [val for val in haploid_s_scores[selection][haploid]["s_minus_median"] if (val >0) ]
            if len(positive_values)<1:
                positive_values = [0]

            for beta_th in range(1,100):
                beta = positive_values[ int(len(positive_values) *  beta_th /100     )]
                haploid_s_scores[selection][haploid]["beta"][beta_th] = beta



    for selection in data3:
        if selection.split("_")[0] != "Control":
            for bait in  data3[selection]:
                for prey in data3[selection][bait]:
                    for bait_bc_ver in data3[selection][bait][prey]:
                        for prey_bc_ver in data3[selection][bait][prey][bait_bc_ver]:
                            for bfg in data3[selection][bait][prey][bait_bc_ver][prey_bc_ver]:
                                ds       = {}
                                ds_med   = data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["s"]+1.0/haploid_s_scores[selection]["Prey_%s_%s"%(prey,prey_bc_ver)]["median"]+1.0/haploid_s_scores[selection]["Bait_%s_%s"%(bait,bait_bc_ver)]["median"]+1.0
                                if method =="PCA":
                                    ds["aMedian"] = ds_med
                                if method =="Y2H":
                                    ds["aMedian"] = ds_med
                                    for b_th in range(1,100):
                                        beta     =  haploid_s_scores[selection]["Bait_%s_%s"%(bait,bait_bc_ver)]["beta"][b_th]
                                        ds["Bth_%03d"%(b_th)] = data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["s"]/beta
                                        if ds["Bth_%03d"%(b_th)] < beta:
                                            ds["Bth_%03d"%(b_th)] = 1.0

                                data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]["ds"] = ds

    data4 = {}

    ######################################################################################################
    ## Changing the data structure
    ## bait_ID -> prey_ID -> category -> selection -> Norm -> 'Raw_data'

    for selection in data3:
        if selection.split("_")[0] != "Control":
            for bait in  data3[selection]:
                for prey in data3[selection][bait]:
                    try:
                        data4[bait][prey][selection] = {}
                    except KeyError:
                        try:
                            data4[bait][prey] = {}
                            data4[bait][prey][selection] = {}
                        except KeyError:
                            data4[bait] = {}
                            data4[bait][prey] = {}
                            data4[bait][prey][selection] = {}

                    for bait_bc_ver in data3[selection][bait][prey]:
                        for prey_bc_ver in data3[selection][bait][prey][bait_bc_ver]:
                            bc = "%s:%s"%(bait_bc_ver,prey_bc_ver)
                            #print bait,bait_bc_ver,prey,prey_bc_ver,bc
                            data4[bait][prey][selection][bc] = {}
                            for bfg in data3[selection][bait][prey][bait_bc_ver][prey_bc_ver]:
                                data4[bait][prey][selection][bc][bfg] = data3[selection][bait][prey][bait_bc_ver][prey_bc_ver][bfg]

    return data4,haploid_s_scores


def output_norm_score(PCA_ds,out_dir):
    scores_LL = {}

    for bait in tqdm(PCA_ds):
        bait_orf =  list(set([i["ORF_SYMBOL"] for i in PCA_db[bait.split("_")[0]] if (i['ID']==bait)] ))[0]
        for prey in PCA_ds[bait]:
            prey_orf = list(set([i["ORF_SYMBOL"] for i in PCA_db[prey.split("_")[0]] if (i['ID']==prey)]))[0]
            for replicate in PCA_ds[bait][prey]:
                ds_array = {}
                for bc_rep in PCA_ds[bait][prey][replicate]:
                    for bfg in PCA_ds[bait][prey][replicate][bc_rep]:
                        for Norm_Method in PCA_ds[bait][prey][replicate][bc_rep][bfg]['ds']:
                            try:
                                ds_array[Norm_Method].append(PCA_ds[bait][prey][replicate][bc_rep][bfg]['ds'][Norm_Method])
                            except KeyError:
                                ds_array[Norm_Method] = [PCA_ds[bait][prey][replicate][bc_rep][bfg]['ds'][Norm_Method] ]
                for Norm_Method in ds_array:
                    try:
                        ds_array[Norm_Method].sort(reverse=True)
                        l = ["PCA",replicate.split("_")[0],replicate.split("_")[1],bait,bait_orf,prey,prey_orf,Norm_Method,np.average(ds_array[Norm_Method]),np.median(ds_array[Norm_Method])]
                        l += ds_array[Norm_Method]
                        #print (l)
                        scores_LL[Norm_Method].append(l)
                    except KeyError:
                        scores_LL[Norm_Method] = [["Method","Selection_condition","Screening_rep","Bait","Bait_ORF","Prey","Prey_ORF","Norm_Method","Average","Median","RANK1","RANK2","RANK3","RANK4","RANK5","RANK6","RANK7","RANK8"]]

    for Norm_Method in scores_LL:
        LL2csv(scores_LL[Norm_Method],"%s/BFG-PCA_scores_Normalized.csv"%(out_dir))

def output_stats(data,db,out_dir):
    #############################
    ## Outping stats as csv in 3_Stats_for_plot
    #############################
    #########################
    ## List of files to output
    # 1. strain_abundance_control.csv   [Method,Condtion,Replicate,bait,bait_ORF,prey,prey_ORF,Raw,F
    # 2. normalized_scores           =  [['Method','Condtion','Replicate','bait','bait_ORF','prey','prey_ORF','s','ds']]
    # 3. autoactivity_level.csv [Method,Strain_type,Conditon,AA_median,AA_median_rank,Strain]
    # 4. bfg_corr.csv           [Method,Conditon,score_type,UpUp,DnDn]
    # 5. diploid_corr.csv       [Method,Conditon,score_type,Diploid1,Diploid2]
    # 6. screening_rep_corr.csv [Method,Conditon,score_type,Rep1,Rep2]
    # 7. ORF_ori_corr.csv       [Method,Conditon,score_type,bait_Prey,Prey_Bait] #Average interanal replicates + screeening replicates
    # 8. screening_meth_corr.csv[score_type,CondA,CondB,ConditonA_score,ConditonB_score] #Average interanal replicates+ screeening replicates
    # beta normalized for each strain, and then average depending on beta -> MCC and choose best beta

    strain_abundance_control    =  [['Method','Replicate','bait','bait_ORF','prey','prey_ORF','bait_bc','prey_bc','BFG','Raw','F']]
    F_heatmap                   =  [['Method','Condition','Replicate','bait','bait_ORF','bait_ID','prey','prey_ORF','prey_ID','bait_bc','prey_bc','Raw','F']]

    scores_heatmap              =  [['Method','Condtion','Replicate','bait','bait_ORF','prey','prey_ORF','s','ds']]

    autoactivity_level          =  [['Method','Strain_type','Conditon','Replicate','Strain','AA_median','AA_median_rank']]
    bfg_corr                    =  [['Method','Conditon','Replicate','score_type','UpUp','DnDn']]
    diploid_corr                =  [['Method','Conditon','Replicate','score_type','Diploid1','Diploid2']]
    screening_rep_corr          =  [['Method','Conditon','score_type','Bait','Prey','Rep1','Rep1_std','Rep2','Rep2_std'] ]
    ORF_ori_corr                =  [['Method','Conditon','score_type','Bait_Prey','Bait_Prey_std','Prey_Bait','Prey_Bait_std'] ]#Average interanal replicates + screeening replicates
    screening_meth_corr         =  [['score_type','CondA','CondB','ConditonA_score','ConditonB_score']] #Average interanal replicates+ screeening replicates

    order = {}
    dhfr12 = [ "%s-%s"%(i["ID"]) for i in db if (i["Type"] =="DHFR12")]
    dhfr12 = list(set(dhfr12))
    dhfr12.sort()
    dhfr3 = [ "%s-%s"%(i["ID"]) for i in db if (i["Type"] =="DHFR3")]
    dhfr3 = list(set(dhfr3))
    dhfr3.sort()


    c = 1
    for i in dhfr12:
        order[i] = c
        c+=1
    c = 1
    for i in dhfr3:
        order[i] = c
        c+=1


    for cond in data:
        F_heatmap  =  [['Method','Condition','Replicate','bait','bait_ORF','bait_ID','prey','prey_ORF','prey_ID','bait_bc','prey_bc','Raw','F']]
        condition = cond.split("_")[0]
        replicate = cond.split("_")[1]
        if condition =="Control":
            for bait in data[cond]:
                bait_ID    =  ("_").join(bait.split("_")[1:])
                bait_order =  order[bait_ID]

                bait_orf =  list(set([i["ORF_SYMBOL"] for i in db[bait.split("_")[0]] if (i['ID']==bait)] ))[0]
                for prey in data[cond][bait]:
                    prey_orf =  list(set([i["ORF_SYMBOL"] for i in db[prey.split("_")[0]] if (i['ID']==prey)] ))[0]
                    prey_ID    =  ("_").join(prey.split("_")[1:])
                    prey_order =  order[prey_ID]

                    for bait_rep in data[cond][bait][prey]:
                        b_order = "%s%s"%(bait_order,bait_rep)
                        for prey_rep in data[cond][bait][prey][bait_rep]:
                            p_order = "%s%s"%(prey_order,prey_rep)
                            ave_raw = []
                            ave_F   = []
                            for bfg in data[cond][bait][prey][bait_rep][prey_rep]:
                                raw = data[cond][bait][prey][bait_rep][prey_rep][bfg]["raw"]
                                F   = data[cond][bait][prey][bait_rep][prey_rep][bfg]["F"]
                                ave_raw.append(raw)
                                ave_F.append(F)


                                l = [method,replicate,bait,bait_orf,prey,prey_orf,bait_rep,prey_rep,bfg,raw,F]
                                strain_abundance_control.append(l)
                            l = [method,replicate,bait,bait_orf,prey,prey_orf,bait_rep,prey_rep,'Ave',np.average(ave_raw),np.average(ave_F)]
                            strain_abundance_control.append(l)

                            l2 = [method,condition,replicate,bait,bait_orf,b_order,prey,prey_orf,p_order,bait_rep,prey_rep,np.average(ave_raw),np.average(ave_F)]
                            F_heatmap.append(l2)
            if not os.path.isdir('%s/heatmap_F'%(out_dir)):
                os.makedirs('%s/heatmap_F'%(out_dir))
            LL2csv(F_heatmap,"%s/F_heatmap_%s_%s.csv"%(out_dir,condition,replicate))

        else:
            for bait in data[cond]:
                bait_orf =  list(set([i["ORF_SYMBOL"] for i in db[bait.split("_")[0]] if (i['ID']==bait)] ))[0]
                bait_ID    =  ("_").join(bait.split("_")[1:])
                bait_order =  order[bait_ID]
                for prey in data[cond][bait]:
                    prey_orf =  list(set([i["ORF_SYMBOL"] for i in db[prey.split("_")[0]] if (i['ID']==prey)] ))[0]
                    prey_ID    =  ("_").join(prey.split("_")[1:])
                    prey_order =  order[prey_ID]

                    for bait_rep in data[cond][bait][prey]:
                        b_order = "%s%s"%(bait_order,bait_rep)
                        for prey_rep in data[cond][bait][prey][bait_rep]:
                            p_order = "%s%s"%(prey_order,prey_rep)
                            ave_raw = []
                            ave_F   = []
                            for bfg in data[cond][bait][prey][bait_rep][prey_rep]:
                                raw = data[cond][bait][prey][bait_rep][prey_rep][bfg]["raw"]
                                F   = data[cond][bait][prey][bait_rep][prey_rep][bfg]["F"]
                                ave_raw.append(raw)
                                ave_F.append(F)
                            l2 = [method,condition,replicate,bait,bait_orf,b_order,prey,prey_orf,p_order,bait_rep,prey_rep,np.average(ave_raw),np.average(ave_F)]
                            F_heatmap.append(l2)
            if not os.path.isdir('%s/heatmap_F'%(out_dir)):
                os.makedirs('%s/heatmap_F'%(out_dir))
            LL2csv(F_heatmap,"%s/F_heatmap_%s_%s.csv"%(out_dir,condition,replicate))

    LL2csv(strain_abundance_control,"%s/strain_abundance_control.csv"%(out_dir))

    for cond in data:
        condition = cond.split("_")[0]
        replicate = cond.split("_")[1]
        D[cond] = {}
        for haploid in data[cond]:
            category = haploid.split("_")[1]
            if category not in D[cond].keys():
                D[cond][category] = []
            #print(data[cond][haploid])
            med = data[cond][haploid]['median']

            l = [method,category,condition,replicate,haploid,med]
            D[cond][category].append(l)
    for cond in D:
        condition = cond.split("_")[0]
        replicate = cond.split("_")[1]
        for category in D[cond]:
            D2 = D[cond][category]
            D3 = sorted(D2,key=lambda x: x[5],reverse=True)
            #pprint(D3)
            i = 1
            for L in D3:
                L2 = L
                L2.append(i)
                autoactivity_level.append(L2)
                #print (L2)
                #print(autoactivity_level[-2:])
                i+=1

    LL2csv(autoactivity_level,"%s/autoactivity_level.csv"%(out_dir))

    diploid_corr  =  [['Method','Conditon','Replicate','Diploid1','Diploid2']]

    REPS[method] = {}
    for cond in data:
        s_heatmap              =  [['Method','Condition','Replicate','bait','prey','s']]
        condition = cond.split("_")[0]
        replicate = cond.split("_")[1]
        if condition !="Control":
            if condition not in REPS[method]:
                REPS[method][condition] = {}
            print (method,condition,replicate)
            for bait in data[cond]:
                bait_orf =  list(set([i["ORF_SYMBOL"] for i in db[bait.split("_")[0]] if (i['ID']==bait)] ))[0]
                bait_ID    =  ("_").join(bait.split("_")[1:])
                bait_order =  order[bait_ID]
                for prey in data[cond][bait]:
                    prey_orf =  list(set([i["ORF_SYMBOL"] for i in db[prey.split("_")[0]] if (i['ID']==prey)] ))[0]
                    prey_ID    =  ("_").join(prey.split("_")[1:])
                    prey_order =  order[prey_ID]
                    pair = "%s-%s"%(bait_orf,prey_orf)


                    dips = {}
                    for bait_rep in data[cond][bait][prey]:
                        b_order = "%s%s"%(bait_order,bait_rep)
                        for prey_rep in data[cond][bait][prey][bait_rep]:
                            p_order = "%s%s"%(prey_order,prey_rep)

                            ave_s = []
                            ave_ds   = []
                            try:
                                bfg_corr.append([method,condition,replicate,'s',data[cond][bait][prey][bait_rep][prey_rep]['UpUp']['s'],data[cond][bait][prey][bait_rep][prey_rep]['DnDn']['s'] ])
                            except KeyError:
                                pass


                            for bfg in data[cond][bait][prey][bait_rep][prey_rep]:
                                s = data[cond][bait][prey][bait_rep][prey_rep][bfg]["s"]

                                ave_s.append(s)

                                try:
                                    REPS[method][condition][pair][replicate].append(s)

                                except KeyError:
                                    try:
                                        REPS[method][condition][pair][replicate] = [s]
                                    except KeyError:

                                            REPS[method][condition][pair] = {}
                                            REPS[method][condition][pair][replicate] = [s]

                            dips["%s:%s"%(bait_rep,prey_rep)]  = {}
                            dips["%s:%s"%(bait_rep,prey_rep)]['s']  =   np.average(ave_s)
                            l2 = [method,condition,replicate,b_order,p_order,np.average(ave_s)]
                            s_heatmap.append(l2)
                    dips_k = list(dips.keys())
                    for bc,bc2 in itertools.combinations(dips_k, 2):
                        criteria = 0
                        if bc.split(":")[0]!= bc2.split(":")[0]:
                            if bc.split(":")[1]!= bc2.split(":")[1]:
                                criteria += 1
                        if criteria >0:
                            l = [method,condition,replicate ,dips[bc]['s'],dips[bc2]['s']]
                            diploid_corr.append(l)
        if not os.path.isdir('%s/heatmap_s'%(out_dir)):
            os.makedirs('%s/heatmap_s'%(out_dir))
        LL2csv(s_heatmap,"%s/heatmap_s/s_heatmap_%s_%s.csv"%(out_dir,method,condition,replicate))

    screening_rep_corr          =  [['Method','Conditon','ORF-pair',"Replicate1","Replicate2",'Rep1','Rep2'] ]

    ORI = {}
    for method in REPS:
        #print (method)
        ORI[method] = {}
        for condition in REPS[method]:

            #print(method,condition)
            for pair in REPS[method][condition]:
                pair_uni = pair.split("-")
                pair_uni.sort()
                pair_uni_str = ("-").join(pair_uni)


                #print (pair,REPS[method][condition][pair].keys())
                for rep1,rep2 in itertools.combinations(list(REPS[method][condition][pair].keys()),2):
                    #print(pair,rep1,rep2)
                    Rep1 = np.average(REPS[method][condition][pair][rep1])
                    Rep2 = np.average(REPS[method][condition][pair][rep2])

                    l = [method,condition,pair,rep1,rep2,Rep1,Rep2]
                    screening_rep_corr.append(l)


                    scrntype1 = "%s%s"%(method,rep1.split("Rep")[0])
                    scrntype2 = "%s%s"%(method,rep2.split("Rep")[0])
                    if scrntype1 == scrntype2:
                        scrntype = "%s %s"%(scrntype1,condition)
                        #print(method,condition,pair,scrntype1,scrntype2)
                        X_Y_ave = np.average([Rep1,Rep2])
                        #Exclude homomeric interations from this analysis
                        if len(set(pair_uni))>1:
                            if pair_uni_str == pair:
                                ori= "X-Y"
                            else:
                                ori= "Y-X"
                            try:
                                ORI[method][scrntype][pair_uni_str][ori] = X_Y_ave
                            except KeyError:
                                try:
                                    ORI[method][scrntype][pair_uni_str] = {}
                                    ORI[method][scrntype][pair_uni_str][ori] = X_Y_ave
                                except KeyError:
                                    ORI[method][scrntype] = {}
                                    ORI[method][scrntype][pair_uni_str] = {}
                                    ORI[method][scrntype][pair_uni_str][ori] = X_Y_ave





    LL2csv(bfg_corr,"%s/BFG_correlation.csv"%(out_dir))
    LL2csv(diploid_corr,"%s/diploid_correlation.csv"%(out_dir))
    LL2csv(screening_rep_corr,"%s/screening_rep_correlation.csv"%(out_dir))
    ORF_ori_corr                =  [['Method','Conditon','X-Y','Y-X'] ]

    conds = []
    for method in ORI:
        for cond in ORI[method]:
            #conds.append(cond)
            for ppi in ORI[method][cond]:
                if len(ORI[method][cond][ppi]) ==2:
                    l = [method,cond,ORI[method][cond][ppi]["X-Y"],ORI[method][cond][ppi]["Y-X"]]
                    ORF_ori_corr.append(l)

    LL2csv(ORF_ori_corr,"%s/ORF_ori_corr.csv"%(out_dir))
