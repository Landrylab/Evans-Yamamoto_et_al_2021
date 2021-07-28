# coding: utf-8
import ast
import os
import sys
import re
from pprint import pprint
import numpy as np
from tqdm import tqdm












#############################
## Defining small functions
#############################
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
def reading(name):
    with open(name, 'r') as F:
        ob_from_file = F.read().replace('\n','')
        mydata = ast.literal_eval(ob_from_file)
    F.close()
    return mydata
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r
def csv2LL(f_name):
    LL = []
    with open(f_name,"r") as F:

        for line in F:
            cols = line.split("\n")[0].split(",")
            LL.append(cols)
    return LL
def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    # This script outputs dict as string, so do not print anything
    print("Generated : %s"%(name))
    F.close()


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

def parse_scores(f,mcc_dir):
    out  = [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","s","Rank","Precision","Recall","MCC",'1_Two-hybrid','2_PCA','3_Union','4_ALL', 'score','TH']]
    DATA = csv2LL(f)
    screen = f.split(".")[0]
    for l in DATA[1:]:
        n      = int(l[3])
        scores = [float(i) for i in l[15:]]
        pair   = l[2]
        sets   = []
        if screen.split("_")[-2] == "Percentile":
            percentile_th = int(screen.split("_")[-1])
            TH =0
            #print(l)
            #if n>32:
            #    #print(n,scores)
            #    for i in range(n):
            #        for j in range(i+1,n):
            #            if scores[i] == scores[j]:
            #                sets.append(scores[i])

            percentile_score = scores[ int(len(scores) *  percentile_th /100     )]
            for score in scores:
                if score == percentile_score:
                    TH = 1
                L = l[:15]
                L.append(score)
                L.append(TH)
                out.append(L)
                TH = 0

        elif screen.split("_")[-1] == "Average":
            Scores = scores
            Scores.append(np.average(scores))
            for score in Scores:
                if score == np.average(scores):
                    TH = 1
                L = l[:15]
                L.append(score)
                L.append(TH)
                out.append(L)
                TH = 0

        elif screen.split("_")[-1] == "Median":
            Scores = scores
            Scores.append(np.median(scores))
            for score in Scores:
                if score == np.median(scores):
                    TH = 1
                L = l[:15]
                L.append(score)
                L.append(TH)
                out.append(L)
                TH = 0



    LL2csv(out,"%s/%s_for_plot.csv"%(mcc_dir,f.split("/")[-1].split(".")[0]))

def biogrid_subset(biogrid,list_of_pairs):
    d = {}
    for pair in list_of_pairs:
        if pair in biogrid:
            d[pair] = biogrid[pair]
    return d

def parse_biogirdtab3(name,taxids):
    D = {"ij":{}}
    binary_methods = ["Two-hybrid",'PCA', 'Biochemical Activity', 'Affinity Capture-Luminescence', 'Reconstituted Complex', 'Co-crystal Structure', 'FRET']

    with open(name,"r") as F:
        for line in tqdm(F):
            #print(line)
            cols = line.split("\n")[0].split("\t")
            #print (cols[15],cols[16])
            """
            0 #BioGRID Interaction ID
            1 Entrez Gene Interactor A
            2 Entrez Gene Interactor B
            3 BioGRID ID Interactor A
            4 BioGRID ID Interactor B
            5 Systematic Name Interactor A
            6 Systematic Name Interactor B
            7 Official Symbol Interactor A
            8 Official Symbol Interactor B
            9 Synonyms Interactor A
            10 Synonyms Interactor B
            11 Experimental System
            12 Experimental System Type
            13 Author
            14 Pubmed ID
            15 Organism Interactor A
            16 Organism Interactor B
            17 Throughput
            18 Score
            19 Modification
            20 Phenotypes
            21 Qualifications
            22 Tags
            23 Source Database
            """
            if cols[12] == "physical":
                #print(cols[15],cols[16],sp.keys())
                if int(cols[15]) in taxids:
                    if int(cols[16]) in taxids:
                        bait_sp = taxids[int(cols[15])]
                        prey_sp = taxids[int(cols[16])]

                        bait = "%s (%s)" %(cols[7],bait_sp)
                        prey =  "%s (%s)" %(cols[8],prey_sp)
                        method = cols[11]
                        #print(method,"\n",line)
                        if method in binary_methods:
                            #print(bait_sp,prey_sp,"\n",line)
                            #print(method)
                            pair = [bait,prey]
                            pair.sort()
                            P = ("_").join(pair)
                            #print(pair)

                            try:
                                D[P][method] +=1
                            except KeyError:
                                try:
                                    D[P][method] = 1
                                except KeyError:
                                    D[P] = {}
                                    D[P][method] = 1
                        else:
                            pair = [bait,prey]
                            pair.sort()
                            P = ("_").join(pair)
                            #print(pair)

                            try:
                                D[P]['other'] +=1
                            except KeyError:
                                try:
                                    D[P]['other'] = 1
                                except KeyError:
                                    D[P] = {}
                                    D[P]['other'] = 1

        F.close()
        return(D)
def cal_MCC(scores,db):
    known_positives = []
    for pair in db:
        for meth in db[pair]:
            if db[pair][meth] > 0:
                known_positives.append(pair)

    known_positives = set(known_positives)
    #pprint(known_positives)
    num_true   = len(known_positives)
    num_sample = float(len(scores))
    num_false  = num_sample - num_true
    MCCmax = -100000000.0
    TP = 0
    FP = 0
    FN = 0
    TN = 0

    for rank in range(1,len(scores)+1):
        score = scores[rank-1]
        """
        score[ 0] ; "Method"
        score[ 1] ; "Selection"
        score[ 2] ; "ORF_pair"
        score[ 3] ; "n"
        score[ 4] ; "Beta"
        score[ 5] ; "Th_method"
        score[ 6] ; "score"
        score[ 7] ; "Rank"
        score[ 8] ; "Precision"
        score[ 9] ; "Recall"
        score[10] ; "MCC"
        score[11] ; '0_Union'
        score[12] ; '1_Two-hybrid'
        score[13] ; '2_PCA'
        score[14] ; '3_Biochemical Activity'
        score[15] ; '4_Co-crystal Structure'
        score[16] ; '5_FRET'
        score[17] ; '6_Reconstituted Complex'
        score[18] ; '7_Affinity Capture-Luminescence'
        score[19:] ; 'Array of scores'
        """
        #print(rank,score[7])
        called_positive = [i[2] for i in scores[:rank]]
        #print(called_positive)
        called_trues = float(len(called_positive))



        TP = float(len([i for i in called_positive if (i in known_positives)]))
        FP = float(rank - TP)
        FN = float(len(known_positives)-TP)
        TN = float(num_sample - (TP + FP + FN))

        try:
            precision = TP / (TP+FP)
        except ZeroDivisionError:
            precision = 0.0
        try:
            recall    = TP / (TP+FN)
        except ZeroDivisionError:
            recall = 0.0

        ## MCC formula
        a = (TP   +   FP)
        b = (TP   +   FN)
        c = (TN   +   FP)
        d = (TN   +   FN)

        if a == 0:
            a = 1
        if b == 0:
            b = 1
        if c == 0:
            c = 1
        if d == 0:
            d = 1


        base = a*b*c*d
        base **= 0.5

        MCCrank  =  (TP   *   TN)
        MCCrank -=  (FP   *   FN)
        MCCrank /=  base




        if MCCrank > MCCmax:
            MCCmax = MCCrank

        scores[rank-1][8]  = precision
        scores[rank-1][9]  = recall
        scores[rank-1][10] = MCCrank



    return scores,MCCmax

def compute_best_MCC(f_nm,out_dir,biogrid,sp):
    MCCmaxLL = [["Method","Selection","Norm","Percentile","MCCmax"]]
    D = {}
    D_ori = {}

    MCC_plot = ["Method","Condition","Norm_method","Percentile"]
    percentiles = [1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99]

    data         = csv2LL(f_nm)
    #print(method_dir,norm_dir,f)
    D = {}
    for l in tqdm(data[1:]):
        #print(l)
        #'Method', 'Selection_condition', 'Screening_rep', 'Bait', 'Bait_ORF', 'Prey', 'Prey_ORF', 'Norm_Method', 'Average', 'Median', 'RANK1', 'RANK2', 'RANK3', 'RANK4', 'RANK5', 'RANK6', 'RANK7', 'RANK8'
        Method             = l[0]
        selection_condition= l[1]
        screening_rep      = l[2]
        bait_strain        = l[3]

        bait_ORF           = "%s (%s)"%(l[4],sp)
        prey_strain        = l[5]
        prey_ORF           = "%s (%s)"%(l[6],sp)
        Norm_method        = l[7]
        #print(l)
        #print(l[10:])
        ds_array           = l[10:]


        pair = [bait_ORF,prey_ORF]
        pair_ID = ("_").join(sorted(pair))

        #print(pair_ID)
        try:
            D[Method][selection_condition][Norm_method][pair_ID] += ds_array
            #print(Method,selection_condition,pair_ID)
        except KeyError:
            try:
                D[Method][selection_condition][Norm_method][pair_ID] = ds_array
            except KeyError:
                try:
                    D[Method][selection_condition][Norm_method] = {}
                    D[Method][selection_condition][Norm_method][pair_ID] = ds_array
                except KeyError:
                    try:
                        D[Method][selection_condition]              = {}
                        D[Method][selection_condition][Norm_method] = {}
                        D[Method][selection_condition][Norm_method][pair_ID] = ds_array
                    except KeyError:
                        D[Method]                                   = {}
                        D[Method][selection_condition]              = {}
                        D[Method][selection_condition][Norm_method] = {}
                        D[Method][selection_condition][Norm_method][pair_ID] = ds_array




    #print("pairs",len(D[Method][selection_condition][Norm_method]))
    #pprint(D)
    d = {}

    for Method in D:
        d[Method] = {}
        for selection_condition in D[Method]:
            d[Method][selection_condition] = {}
            for Norm_method in  D[Method][selection_condition]:

                #print(Method,selection_condition,Norm_method,len(D[Method][selection_condition][Norm_method]))
                #@pprint(D[Method][selection_condition][Norm_method])
                #print(jso)
                for pair in D[Method][selection_condition][Norm_method]:
                    DS = [ float(i) for i in D[Method][selection_condition][Norm_method][pair]]
                    DS2 = sorted(DS)



                    try:
                        twohyb   = biogrid[pair]['Two-hybrid']
                    except KeyError:
                        twohyb = 0
                    try:
                        pca      = biogrid[pair]['PCA']
                    except KeyError:
                        pca = 0
                    try:

                        biochem  = biogrid[pair]['Biochemical Activity']
                    except KeyError:
                        biochem = 0
                    try:

                        cocrystal = biogrid[pair]['Co-crystal Structure']
                    except KeyError:
                        cocrystal = 0
                    try:

                        FRET     = biogrid[pair]['FRET']
                    except KeyError:
                        FRET = 0
                    try:

                        recomp   = biogrid[pair]['Reconstituted Complex']
                    except KeyError:
                        recomp = 0
                    try:

                        AffWesturn= biogrid[pair]['Affinity Capture-Luminescence']
                    except KeyError:
                        AffWesturn = 0
                    try:

                        other= biogrid[pair]['Affinity Capture-Luminescence']
                    except KeyError:
                        other = 0

                    union    = twohyb + pca+biochem+cocrystal+FRET+recomp+AffWesturn
                    ALL      = union+other
                    #print(pair)
                    #print(DS2)
                    for percentile in percentiles:
                        th_meth = "Percentile_%02d"%(percentile)
                        percentile_score = DS2[ int(len(DS2) *  percentile /100     )]
                        l = [Method,selection_condition,pair,len(DS2),Norm_method,"Percentile_%02d"%(percentile),percentile_score,0,"X","X","X", twohyb ,pca ,union ,ALL]
                        l += DS2
                        #print(th_meth,l)

                        try:
                            d[Method][selection_condition][Norm_method][th_meth].append(l)
                        except KeyError:
                            try:

                                d[Method][selection_condition][Norm_method][th_meth] =[["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'0_Union','1_Two-hybrid','2_PCA', '3_Biochemical Activity','4_Co-crystal Structure', '5_FRET', '6_Reconstituted Complex','7_Affinity Capture-Luminescence']]
                                d[Method][selection_condition][Norm_method][th_meth].append(l)
                            except KeyError:
                                d[Method][selection_condition][Norm_method] = {}
                                d[Method][selection_condition][Norm_method][th_meth] = [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'0_Union','1_Two-hybrid','2_PCA', '3_Biochemical Activity','4_Co-crystal Structure', '5_FRET', '6_Reconstituted Complex','7_Affinity Capture-Luminescence']]
                                d[Method][selection_condition][Norm_method][th_meth].append(l)
                    #print(percentile,len(d[Method][selection_condition][Norm_method][th_meth]))
                    #print(th_meth,d[Method][selection_condition][Norm_method][th_meth])
                    #print(header)
                    #pprint(hoge)

                    th_meth = "Average"
                    l = [Method,selection_condition,pair,len(DS2),Norm_method,"Average",np.average(DS2),0,"X","X","X", twohyb ,pca ,union ,ALL]
                    l += DS2
                    try:
                        d[Method][selection_condition][Norm_method][th_meth].append(l)
                    except KeyError:
                        try:
                            d[Method][selection_condition][Norm_method][th_meth] = [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'0_Union','1_Two-hybrid','2_PCA', '3_Biochemical Activity','4_Co-crystal Structure', '5_FRET', '6_Reconstituted Complex','7_Affinity Capture-Luminescence']]
                            d[Method][selection_condition][Norm_method][th_meth].append(l)
                        except KeyError:
                            d[Method][selection_condition][Norm_method] = {}
                            d[Method][selection_condition][Norm_method][th_meth] = [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'0_Union','1_Two-hybrid','2_PCA', '3_Biochemical Activity','4_Co-crystal Structure', '5_FRET', '6_Reconstituted Complex','7_Affinity Capture-Luminescence']]
                            d[Method][selection_condition][Norm_method][th_meth].append(l)

                    th_meth = "Median"
                    l = [Method,selection_condition,pair,len(DS2),Norm_method,"Median",np.median(DS2),0,"X","X","X", twohyb ,pca ,union ,ALL]
                    l += DS2
                    try:
                        d[Method][selection_condition][Norm_method][th_meth].append(l)
                    except KeyError:
                        try:
                            d[Method][selection_condition][Norm_method][th_meth] [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'1_Two-hybrid','2_PCA', '3_Union','4_ALL']]
                            d[Method][selection_condition][Norm_method][th_meth].append(l)
                        except KeyError:

                            d[Method][selection_condition][Norm_method] = {}
                            d[Method][selection_condition][Norm_method][th_meth] = [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'1_Two-hybrid','2_PCA', '3_Union','4_ALL']]
                            d[Method][selection_condition][Norm_method][th_meth].append(l)
    BEST_MCC = 0
    for method in d:
        for sel in d[method]:
            #print(sel)
            for norm in  d[method][sel]:
                print("MCC      ",method,sel,norm)
                for th_meth in tqdm(d[method][sel][norm]):
                    #print(method,sel,norm,th_meth,len(d[method][sel][norm][th_meth]))

                    ps = [i[2] for i in d[method][sel][norm][th_meth]]
                    biogrid = biogrid_subset(biogrid,ps)

                    '''
                    ########
                    ## Taking 'relieable' (detected in more than tw papers) interactions
                    reliable_dataset = {}
                    for p in biogrid:
                        n = 0
                        for m in biogrid[p]:
                            n += biogrid[p][m]
                        if n  >1:
                            reliable_dataset[p] = biogrid[p]

                    #print(len(set(ps)),len(biogrid))
                    #pprint(biogrid)
                    '''


                    a  = []
                    LL = d[method][sel][norm][th_meth]
                    ##### Add rank info to the LL
                    header = LL[0]
                    scores = LL[1:]

                    scores.sort(key=lambda x: x[6],reverse=True)
                    Rank = 1
                    for score in scores:
                        score[7] = Rank
                        Rank+=1
                    #
                    #pprint(scores)
                    scores,MCCmax = cal_MCC(scores,biogrid)
                    info = [method,sel,norm,th_meth]
                    out = [["Method","Selection","ORF_pair","n","Norm_meth","Th_method","score","Rank","Precision","Recall","MCC",'1_Two-hybrid','2_PCA', '3_Union','4_ALL']]
                    out+= scores
                    out_f = "%s/all_conditions/%s/SCORES_BFG-%s.csv"%(out_dir,sp,("_").join(info))
                    LL2csv(out,out_f)
                    info.append(MCCmax)
                    MCCmaxLL.append(info)
                    #print(info)

                    if MCCmax > BEST_MCC:
                        BEST_performing = out_f
                        BEST_MCC = MCCmax
        LL2csv(MCCmaxLL,"%s/MCCmax_%s.csv"%(out_dir,sp))
    return BEST_performing

def biogrid_subset(biogrid,list_of_pairs):
    d = {}
    for pair in list_of_pairs:
        if pair in biogrid:
            d[pair] = biogrid[pair]
    return d
