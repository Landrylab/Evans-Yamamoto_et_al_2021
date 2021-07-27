import numpy as np
import math
import os

from random import seed
from random import random

PI    = 4 * math.atan2(1, 1)

def randn(m,sigma):
  r1  = random()
  r2  = random()
  while r1==0:
      r1 = random()
  value = (sigma* (( -2 * math.log(r1) )**1/2) * math.sin(2 * PI * r2)) + m
  return value


def lograndn(CV):
  sigma = math.log((CV**2)+1)
  sigma **= 1/2

  return math.exp(randn(0,sigma))


def mean_sigma(array):
    mean = np.mean(array)
    sigma = 0

    for n in array:
        sigma += (n-mean)**2;

    sigma /= len(array)
    sigma **= 1/2

    CV = sigma/mean

    return (mean,sigma,CV)



def adjust_num(array,CV,pop_size):
    for n in range(len(array)):
        val = array[n][1]* lograndn(CV)
        array[n][1] = val

    sum = np.sum([i[1] for i in array])

    for n in range(len(array)):
        val = float(array[n][1])*(float(pop_size)/sum)
        array[n][1] = val
    return array



def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    # This script outputs dict as string, so do not print anything
    print("Generated : %s"%(name))
    F.close()

def hap_coverage(x,y,nm):
    all = [i[1] for i in x]
    all += [i[1] for i in y]
    mx = int(math.log(max(all),10)+2)
    cov = [["Threashold".ljust(mx+1),"x_percentage".ljust(mx+1),"y_percentage".ljust(mx+1)]]
    for th in range(0,mx):
        val_x = above_threshold([i[1] for i in x],10**th)*100
        val_y = above_threshold([i[1] for i in y],10**th)*100
        cov.append([str(10**th).ljust(mx+1),str(round(val_x,2)).ljust(mx+1),str(round(val_y,2)).ljust(mx+1)])
    return cov

def dip_coverage(dip,nm):
    all = [i[1] for i in dip]
    mx = int(math.log(max(all),10)+2)
    cov = [["Threashold".ljust(mx+1),"Percentage".ljust(mx+1)]]
    for th in range(0,mx):
        val = above_threshold([i[1] for i in dip],10**th)*100
        cov.append([str(10**th).ljust(mx+1),str(round(val,2)).ljust(mx+1)])
    return cov

def above_threshold(pop,th):
    above = [ strain for strain in pop if (strain > th)]
    value = float(len(above))/len(pop)
    return value

def print_LL(LL):
    for L in LL:
        l = [str(i) for i in L]
        print(("\t").join(l))





def mating(hap_x,hap_y):
    dip = []
    for strain_x in hap_x:
        for strain_y in hap_y:
            xy_value = strain_x[1] * strain_y[1]
            dip.append([["%s"%(strain_x[0]),"%s"%(strain_y[0])],xy_value])
    return dip

def positive_interaction(array,CV_pos,CV_autoactivity,pos_rate,pop_size):

    strains = [i[0][0] for i in array]
    strains += [i[0][1] for i in array]

    autoactivity = {}
    for strain in strains:
        autoactivity[strain] = lograndn(CV_autoactivity)

    for i in array:
        signal = lograndn(CV_pos)
        val    = i[1] * signal
        i[1]   = val

    array.sort(key = lambda x:x[1],reverse=True)

    positives = len(array)*pos_rate
    for i in range(len(array)):
        if i<positives:
            array[i][0].append(1)
        else:
            array[i][0].append(0)

    for i in array:
        x = i[0][0]
        y = i[0][1]

        aa = 0
        aa = autoactivity[x]*autoactivity[y]

        value = i[1] * aa
        i[1] = value

    sum = np.sum([i[1] for i in array])

    for n in range(len(array)):
        val = float(array[n][1])*(float(pop_size)/sum)
        array[n][1] = val

    return array

def abundance(array):
    dip = {}
    sum = 0
    for i in array:
        dip_nm = ("-").join([i[0][0],i[0][1]])
        try:
            dip[dip_nm]  += i[1]+1
        except KeyError:
            dip[dip_nm]  = i[1]+1
        sum += i[1]+1

    for i in dip:
        dip[i] /= sum

    return dip



def marginal(array):
    x_hap = {}
    y_hap = {}
    sum = 0

    for i in array:
        x = i[0][0]
        y = i[0][1]
        try:
            x_hap[x]  += i[1]+1
        except KeyError:
            x_hap[x]  = i[1]+1
        try:
            y_hap[y]  += i[1]+1
        except KeyError:
            y_hap[y]  = i[1]+1
        sum += i[1]+1

    for i in x_hap:
        x_hap[i] /= sum
    for i in y_hap:
        y_hap[i] /= sum

    return [x_hap[i] for i in x_hap],[y_hap[i] for i in y_hap]

def format_array(array):
    LL = [["x","y","value"]]
    for i in array:
        l = [i[0][0],i[0][1],i[1]]
        LL.append(l)
    return LL
