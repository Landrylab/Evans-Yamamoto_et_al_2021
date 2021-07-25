import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
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
        val = array[n]* lograndn(CV)
        array[n] = val

    sum = np.sum(array)

    for n in range(len(array)):
        val = array[n]/sum*pop_size
        array[n] = val
    return array

def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    # This script outputs dict as string, so do not print anything
    print("Generated : %s"%(name))
    F.close()

def hap_coverage(x,y,nm):
    all = x
    all += y
    mx = int(math.log(max(all),10)+2)
    cov = [["Cells".ljust(mx+1),"x_percentage".ljust(mx+1),"y_percentage".ljust(mx+1)]]
    for th in range(0,mx):
        val_x = above_threshold(x,10**th)*100
        val_y = above_threshold(y,10**th)*100
        cov.append([str(10**th).ljust(mx+1),str(round(val_x,2)).ljust(mx+1),str(round(val_y,2)).ljust(mx+1)])
    return cov


def above_threshold(pop,th):
    above = [ strain for strain in pop if (strain > th)]
    value = float(len(above))/len(pop)
    return value

def print_LL(LL):
    for L in LL:
        l = [str(i) for i in L]
        print(("\t").join(l))
