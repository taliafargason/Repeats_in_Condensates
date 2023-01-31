from Bio import SeqIO
import re
def mean(list):
    try:
        return sum(list)/len(list)
    except ZeroDivisionError:
        return 0
inputfile='reviewed_only_organism9606_existence1_fragmentfalse_reviewedtrue_2022.08.17.fasta'
seq=list(SeqIO.parse(inputfile, "fasta"))
phase_list=[]
rrmlist = open("proteins_w_rrm_domains.txt","r")
RRMlist = [re.split('\s+',line)[0] for line in rrmlist]
with open ('combo.txt') as p:
    for line in p:
        if line[:-1] not in phase_list:
            phase_list.append(line[:-1])
 
 
amino_acids = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']

#create a list of all possible dipeptide combinations
dipeptides = []
for x in range(0,len(amino_acids)):
        a = amino_acids[x]
        for n in range(x,(len(amino_acids))):
            b = amino_acids[n]
            dipeptides.append([[a],[b]])


#create a dictionary that has a list of all possible combinations for each peptide
lengths = {}
fract_in_condensates = {}
fract_in_condensates_w_RRM = {}
num = {}
num_w_RRM = {}
avg_in_condensates = {}
avg_in_condensates_w_RRM = {}
avg_not_in_condensates = {}
avg_not_in_condensates_w_RRM = {}

for y in dipeptides:
    for a in y[0]:
        for b in y[1]:
            lengths['%s%s'%(a,b)] = []
            fract_in_condensates['%s%s'%(a,b)] = []
            fract_in_condensates_w_RRM['%s%s'%(a,b)] = []
            num['%s%s'%(a,b)] = []
            num_w_RRM['%s%s'%(a,b)] = []
            avg_in_condensates['%s%s'%(a,b)] = []
            avg_in_condensates_w_RRM['%s%s'%(a,b)] = []
            avg_not_in_condensates['%s%s'%(a,b)] = []
            avg_not_in_condensates_w_RRM['%s%s'%(a,b)] = []
            c =(b+a)
            d =(a+b)
            search_motif1 = [c,d]

            #construct a list of proteins that don't have 2-mer repeats
            namelist1 = []
            for i in range(len(seq)):
                namelist2 = []
                for j in range(len(search_motif1)):
                    if (search_motif1[j] in seq[i].seq and seq[i].id.split('|')[1] not in namelist2):
                        namelist2.append(seq[i].id.split('|')[1])
                if seq[i].id.split('|')[1] not in namelist2:
                    namelist1.append(seq[i].id.split('|')[1])
            proteins_w_repeat = 0
            proteins_w_repeat_in_condensates = 0
            proteins_w_repeat_n_RRM = 0
            proteins_w_repeat_n_RRM_in_condensates = 0
            #The lists below only apply to proteins that contain the repeat type:
            #namelist = proteins that contain the repeat
            judge=[]
            #judge = whether the protein is in condensates
            RRM_RS_list = []
            #RRM_RS_list = whether the protein contains an RRM domain    
            for v in range(len(namelist1)):
                proteins_w_repeat +=1
                if namelist1[v] in phase_list:
                    proteins_w_repeat_in_condensates +=1
                    judge.append(1)
                    if namelist1[v] in RRMlist:
                        RRM_RS_list.append(1)
                        proteins_w_repeat_n_RRM +=1
                        proteins_w_repeat_n_RRM_in_condensates +=1
                    else:
                        RRM_RS_list.append(0)
                else:
                    judge.append(0)
                    if namelist1[v] in RRMlist:
                        RRM_RS_list.append(1)
                        proteins_w_repeat_n_RRM +=1
                    else:
                        RRM_RS_list.append(0)
            
            try:
                fraction_in_condensates_total = proteins_w_repeat_in_condensates/proteins_w_repeat
            except ZeroDivisionError:
                fraction_in_condensates_total = 'N/A'
            
            fract_in_condensates['%s%s'%(a,b)] += ["%s"%fraction_in_condensates_total]
            num['%s%s'%(a,b)] += ["%s"%proteins_w_repeat]
            try:
                fraction_in_condensates_w_RRM = proteins_w_repeat_n_RRM_in_condensates/proteins_w_repeat_n_RRM
            except ZeroDivisionError:
                fraction_in_condensates_w_RRM = 'N/A'
            
            fract_in_condensates_w_RRM['%s%s'%(a,b)] += ["%s"%fraction_in_condensates_w_RRM]
            num_w_RRM['%s%s'%(a,b)] += ["%s"%proteins_w_repeat_n_RRM]
            
        #construct lists of proteins that have 2-mer to 8-mer repeats
            for u in range(1,5):
                c =((b+a)*u)
                d =((a+b)*u)
                #e = [c,d]
                lengths['%s%s'%(a,b)] += [[c,d]]
                search_motif = [c,d]
                #lengths.append(e)  
                #we want to print side-by-side 8-membered lists of % in condensates:
                #calculate % in condensates for each item in the dictionary
                #the input for the function needs to be fixed in order to get the correct p-values.  Totals are not carried over.  May need to do the first iteration separately so totals can be inserted.
                #g = open('RRM_RS_report_%s.txt'%search_motif[0],'w')
                #seq=list(SeqIO.parse(inputfile, "fasta"))
                #list of proteins containing RS8 and SR8
                name_list=[]
                count_list=[]
                count_list_ps =[]
                count_list_RRM =[]
                count_list_ps_RRM =[]
                assessment1 = []
                assessment2 = []
                for i in range(len(seq)):
                    for j in range(len(search_motif)):
                        if (search_motif[j] in seq[i].seq) and (seq[i].id.split('|')[1] not in name_list):
                            #if the protein contains repeat and is not already in namelist, append uniprotID to namelist
                            name_list.append(seq[i].id.split('|')[1])

                            #determine the number of times a repeat of the given length appears
                            temp=[]
                            for k in range(len(search_motif)):
                                temp.append(seq[i].seq.count(search_motif[k]))
                            number_of_repeats=max(temp)
                            count_list.append(number_of_repeats)
                            if seq[i].id.split('|')[1] in RRMlist:
                                count_list_RRM.append(number_of_repeats)
                            if seq[i].id.split('|')[1] in phase_list:
                                count_list_ps.append(number_of_repeats)
                                if seq[i].id.split('|')[1] in RRMlist:
                                    count_list_ps_RRM.append(number_of_repeats)
                            #number times search motif appears


                avg_count = mean(count_list)
                avg_count_RRM = mean(count_list_RRM)
                avg_count_ps = mean(count_list_ps)
                avg_count_ps_RRM = mean(count_list_ps_RRM)              
                avg_in_condensates['%s%s'%(a,b)] += ["%s"%avg_count_ps]
                avg_in_condensates_w_RRM['%s%s'%(a,b)] += ["%s"%avg_count_ps_RRM]
                avg_not_in_condensates['%s%s'%(a,b)] += ["%s"%avg_count]
                avg_not_in_condensates_w_RRM['%s%s'%(a,b)] += ["%s"%avg_count_RRM]
                
                proteins_w_repeat = 0
                proteins_w_repeat_in_condensates = 0
                proteins_w_repeat_n_RRM = 0
                proteins_w_repeat_n_RRM_in_condensates = 0
                #The lists below only apply to proteins that contain the repeat type:
                #namelist = proteins that contain the repeat
                judge=[]
                #judge = whether the protein is in condensates
                RRM_RS_list = []
                #RRM_RS_list = whether the protein contains an RRM domain
                
                for v in range(len(name_list)):
                    proteins_w_repeat +=1
                    if name_list[v] in phase_list:
                        proteins_w_repeat_in_condensates +=1
                        judge.append(1)
                        if name_list[v] in RRMlist:
                            RRM_RS_list.append(1)
                            proteins_w_repeat_n_RRM +=1
                            proteins_w_repeat_n_RRM_in_condensates +=1
                        else:
                            RRM_RS_list.append(0)
                    else:
                        judge.append(0)
                        if name_list[v] in RRMlist:
                            RRM_RS_list.append(1)
                            proteins_w_repeat_n_RRM +=1
                        else:
                            RRM_RS_list.append(0)
                
                try:
                    fraction_in_condensates_total = proteins_w_repeat_in_condensates/proteins_w_repeat
                except ZeroDivisionError:
                    fraction_in_condensates_total = 'N/A'
                
                fract_in_condensates['%s%s'%(a,b)] += ["%s"%fraction_in_condensates_total]
                num['%s%s'%(a,b)] += ["%s"%proteins_w_repeat]
                try:
                    fraction_in_condensates_w_RRM = proteins_w_repeat_n_RRM_in_condensates/proteins_w_repeat_n_RRM
                except ZeroDivisionError:
                    fraction_in_condensates_w_RRM = 'N/A'
                
                fract_in_condensates_w_RRM['%s%s'%(a,b)] += ["%s"%fraction_in_condensates_w_RRM]
                num_w_RRM['%s%s'%(a,b)] += ["%s"%proteins_w_repeat_n_RRM]
                #average_number_repeats
                #average_number_repeats_in_condensates
                #average_number_repeats_w_RRM
                #average_number_repeats_w_RRM_in_condensates
                            
total_fract = open('fraction_in_condensates_total_d.txt','w')
total_fract.write('no_2mer 2mer 4mer 6mer 8mer\n')
total = open('total_d.txt','w')
total.write('no_2mer 2mer 4mer 6mer 8mer\n')
RRM_fract = open('fraction_proteins_with_RRM_in_condensates_total_d.txt','w')
RRM_fract.write('no_2mer 2mer 4mer 6mer 8mer\n')
total_RRM = open('total_RRM_d.txt','w')
total_RRM.write('no_2mer 2mer 4mer 6mer 8mer\n')
for line in fract_in_condensates:
         total_fract.write(line)
         print(fract_in_condensates[line],file = total_fract)
         #total_fract.write(" %s %s %s %s\n"%(fract_in_condensates[line][0],fract_in_condensates[line][1],fract_in_condensates[line][2],fract_in_condensates[line][3]))
for line in fract_in_condensates_w_RRM:
         RRM_fract.write(line)
         print(fract_in_condensates[line],file = RRM_fract)
for line in num:
    total.write(line)
    print(num[line],file=total)
for line in num_w_RRM:
    total_RRM.write(line)
    print(num_w_RRM[line],file=total_RRM)
         #RRM_fract.write(" %s %s %s %s\n"%(fract_in_condensates_w_RRM[line][0],fract_in_condensates_w_RRM[line][1],fract_in_condensates_w_RRM[line][2],fract_in_condensates_w_RRM[line][3]))
#overall.write("in_condensates not_in_condensates\n")  
total.close()
total_RRM.close()
total_fract.close()
RRM_fract.close()

ave = open('ave_total_d.txt','w')
ave.write('2mer 4mer 6mer 8mer\n')
RRM_ave = open('ave_proteins_with_RRM_total_d.txt','w')
RRM_ave.write('2mer 4mer 6mer 8mer\n')
ave_ps = open('ave_in_condensates_total_d.txt','w')
ave_ps.write('2mer 4mer 6mer 8mer\n')
RRM_ave_ps = open('ave_proteins_with_RRM_in_condensates_d.txt','w')
RRM_ave_ps.write('2mer 4mer 6mer 8mer\n')

for line in avg_not_in_condensates:
         ave.write(line)
         print(avg_not_in_condensates[line],file = ave)
         #ave.write(" %s %s %s %s\n"%(avg_not_in_condensates[line][0],avg_not_in_condensates[line][1],avg_not_in_condensates[line][2],avg_not_in_condensates[line][3]))
for line in avg_not_in_condensates_w_RRM:
         RRM_ave.write(line)
         print(avg_not_in_condensates_w_RRM[line],file = RRM_ave)
         #RRM_ave.write(" %s %s %s %s\n"%(avg_not_in_condensates_w_RRM[line][0],avg_not_in_condensates_w_RRM[line][1],avg_not_in_condensates_w_RRM[line][2],avg_not_in_condensates_w_RRM[line][3]))
for line in avg_in_condensates:
         ave_ps.write(line)
         print(avg_in_condensates[line], file = ave_ps)
         #ave_ps.write(" %s %s %s %s\n"%(avg_in_condensates[line][0],avg_in_condensates[line][1],avg_in_condensates[line][2],avg_in_condensates[line][3]))
for line in avg_in_condensates_w_RRM:
         RRM_ave_ps.write(line)
         print(avg_in_condensates_w_RRM[line],file = RRM_ave)
         #RRM_ave_ps.write(" %s %s %s %s\n"%(avg_in_condensates_w_RRM[line][0],avg_in_condensates_w_RRM[line][1],avg_condensates_w_RRM[line][2],avg_in_condensates_w_RRM[line][3]))

#overall.write("in_condensates not_in_condensates\n")  
ave.close()
RRM_ave.close()
ave_ps.close()
RRM_ave_ps.close()
total_fract.close()
RRM_fract.close()

