# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 18:05:00 2022

@author: tfarg
"""
size = 6
repeat_max = 5
import re
import random
import numpy as np
from scipy.stats import fisher_exact
#open the file for final percentages
smp = open('1-13-23_size_matched_percentages.txt','w')


lengths = []
#create a list of all our desired dipeptides at all possible lengths
for a in 'R':
    for b in 'S':
        for i in range(1,repeat_max):
            c =((b+a)*i)
            d =((a+b)*i)
            e = [c,d]
            lengths.append(e)
                    
#open each file in the list and create a new outfile simultaneously
for search_motif in lengths:
    g = open('RRM_RS_report_%s.txt'%search_motif[0],'r')
    sml = open('1-13-23-size_matched_lists_%s.txt'%search_motif[0],"w")
    proteins_w_repeat = [re.split("\s+",line) for line in g]
    g.close()
    smp.write("\n%s_RRM no_%s_RRM %s_no_RRM no_%s_no_RRM\n"%(search_motif[0],search_motif[0],search_motif[0],search_motif[0]))
    #for each new file containing the protein info, separate proteins w/RRM from w/o RRM and w/ repeat from w/o repeat
    repeat_no_RRM = []
    repeat_RRM = []
    no_repeat_RRM = []
    no_repeat_no_RRM = []
    for line in proteins_w_repeat:
            if len(line)==5:
                ID = line[0]
                repeat = int(line[1])
                phase = int(line[2])
                RRM = int(line[3])
                if repeat > 0 and RRM > 0:
                    repeat_RRM.append([ID,phase])
                elif repeat > 0 and RRM == 0:
                    repeat_no_RRM.append([ID,phase])
                elif repeat == 0 and RRM > 0:
                    no_repeat_RRM.append([ID,phase])
                elif repeat == 0 and RRM == 0:
                    no_repeat_no_RRM.append([ID,phase])
                else: 
                    print ("ERROR %s"%line)
            else: 
                print ("ERROR %s"%line)
                
    #now that we have read in our files and sorted into four lists, create new lists with 50 sets of random size-matched subsets
    rand_rep_RRM = []
    rand_norep_RRM = []
    rand_rep_noRRM = []
    rand_norep_noRRM = []
    
    for n in range(0,50):
        try:
            repeat_RRM_size_match = [random.sample(repeat_RRM,size)]
            repeat_RRM_0 = 0
            repeat_RRM_nops = 0
        except ValueError:
            repeat_RRM_size_match = [random.sample(repeat_RRM,len(repeat_RRM))]
            repeat_RRM_0 = 0
            repeat_RRM_nops = 0

        try:
            repeat_no_RRM_size_match = [random.sample(repeat_no_RRM,size)]
            repeat_no_RRM_0 = 0
            repeat_no_RRM_nops = 0
        except ValueError:
            repeat_no_RRM_size_match = [random.sample(repeat_no_RRM,len(repeat_no_RRM))]
            repeat_no_RRM_0 = 0
            repeat_no_RRM_nops = 0
        try:
            no_repeat_RRM_size_match = [random.sample(no_repeat_RRM,size)]
            no_repeat_RRM_0 = 0
            no_repeat_RRM_nops = 0
        except ValueError:
            no_repeat_RRM_size_match = [random.sample(no_repeat_RRM,len(no_repeat_RRM))]
            no_repeat_RRM_0 = 0
            no_repeat_RRM_nops = 0
        try:
            no_repeat_no_RRM_size_match = [random.sample(no_repeat_no_RRM,size)]
            no_repeat_no_RRM_0 = 0
            no_repeat_no_RRM_nops = 0
        except ValueError:
            no_repeat_no_RRM_size_match = [random.sample(no_repeat_no_RRM,len(no_repeat_no_RRM))]
            no_repeat_no_RRM_0 = 0
            no_repeat_no_RRM_nops = 0
        sml.write("\n%s_RRM no_%s_RRM %s_no_RRM no_%s_no_RRM\n"%(search_motif[0],search_motif[0],search_motif[0],search_motif[0]))
        for o in range(0,size):
            sml.write("%s "%(repeat_RRM_size_match[0][o][0]))
            if repeat_RRM_size_match[0][o][1] ==1:
                repeat_RRM_0 +=1
            else:
                repeat_RRM_nops +=1
            sml.write("%s "%(no_repeat_RRM_size_match[0][o][0]))
            if no_repeat_RRM_size_match[0][o][1] ==1:
                no_repeat_RRM_0 +=1
            else:
                no_repeat_RRM_nops +=1
            sml.write("%s "%(repeat_no_RRM_size_match[0][o][0]))
            if repeat_no_RRM_size_match[0][o][1] ==1:
                repeat_no_RRM_0 +=1
            else:
                repeat_no_RRM_nops +=1
            sml.write("%s\n"%(no_repeat_no_RRM_size_match[0][o][0]))
            if no_repeat_no_RRM_size_match[0][o][1] ==1:
                no_repeat_no_RRM_0 +=1
            else:
                no_repeat_no_RRM_nops +=1
        try:
            repeat_RRM_percentage = int(repeat_RRM_0)/int(size)*100
            smp.write("%s "%repeat_RRM_percentage)
        except ZeroDivisionError:
            smp.write(" ")
        try:
            norepeat_RRM_percentage = int(no_repeat_RRM_0)/int(size)*100
            smp.write("%s "%norepeat_RRM_percentage)
        except ZeroDivisionError:
            smp.write(" ")
        try:
            repeat_noRRM_percentage = int(repeat_no_RRM_0)/int(size)*100
            smp.write("%s "%repeat_noRRM_percentage)
        except ZeroDivisionError:
            smp.write(" ")
        try:    
            norepeat_noRRM_percentage = int(no_repeat_no_RRM_0)/int(size)*100
            smp.write("%s \n"%norepeat_noRRM_percentage)
        except ZeroDivisionError:
            smp.write(" \n")
        sml.write("\n\n")
    sml.close()
smp.close()

                    
                
