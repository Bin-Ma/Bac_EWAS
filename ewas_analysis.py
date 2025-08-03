#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 18:52:13 2025

@author: mabin
"""


###EWAS analysis 

import pandas  as pd 
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('index',help="Strain sequence index location.")
parser.add_argument('gene_index',help="Gene index location.")
parser.add_argument('out',help="Output location.")
parser.add_argument('methy_info',help="Pan-Methylome location.")
parser.add_argument('gene_info',help="Gene information (Gene_name,Location,Length).")
parser.add_argument('classA',help="Strain name adapted to the environment.")


parser.add_argument('--methy_name', default='6mA' , help="Methylation name (6mA of 5mC), default = 6mA.")
parser.add_argument('--methy_type', default='complete' , help="Methylation type (complete or pos or neg), default = complete.")
parser.add_argument('--gene_type', default='G' , help="Gene type: gene(G) or non-coding(B), default = G.")

args = parser.parse_args()


##input ess gene

seq_index = list(pd.read_csv(args.gene_index,header = None)[0])
for i in range(len(seq_index)):
    seq_index[i] =  seq_index[i].split('.fa')[0]

seq_index = args.index

NGS_info = list(pd.read_csv(seq_index,header = None)[0])

methy_info = pd.read_csv(args.methy_info,header=None)

core_gene_info = pd.read_csv(args.gene_info)



##gene site NGS_info
methy_match = pd.DataFrame({'Gene':[],'site':[]})
methy_match_NGS = []


def methy_info_match(methy_info,gene_info,NGS_info,methy_type='6mA',methy_loc='complete',gene_type='G'):
    
    for i in range(len(methy_info)):
        if pd.isna(methy_info[7][i]):
            continue
        
        gene_name = gene_type+'_'+methy_info[0][i]
        if gene_name in list(gene_info['Gene']):
            gene_index = list(gene_info['Gene']).index(gene_name)
        else:
            continue
        
        gene_length = gene_info['length'][gene_index]
        gene_start = gene_info['start'][gene_index]

        if methy_loc == 'complete':
            ##pos
            methy_site = methy_info[1][i]
            if len(methy_match[(methy_match['Gene']==gene_name) & (methy_match['site']==methy_site)]) == 0:
                methy_match.loc[len(methy_match)] =  {'Gene':gene_name,'site':methy_site}
                tmp_methy = []
                for j in range(len(NGS_info)):
                    if NGS_info[j] in methy_info[7][i].split(','):
                        if methy_type == '6mA':
                            #tmp_methy.append('GM')
                            tmp_methy.append('AMMA')
                        else:
                            #tmp_methy.append('CM')
                            tmp_methy.append('CNNC')

                    else:
                        if methy_type == '6mA':
                            tmp_methy.append('AAAA')
                        else:
                            tmp_methy.append('CCCC')
                        
                methy_match_NGS.append(tmp_methy)
                
            else:
                methy_index = methy_match[(methy_match['Gene']==gene_name) & (methy_match['site']==methy_site)].index[0]
                for j in range(len(NGS_info)):
                    if methy_match_NGS[methy_index][j] == 'AAAA'  or  methy_match_NGS[methy_index][j] == 'AMAA'  or methy_match_NGS[methy_index][j] == 'AAMA'  or methy_match_NGS[methy_index][j] == 'CCCC'  or methy_match_NGS[methy_index][j] == 'CNCC'  or methy_match_NGS[methy_index][j] == 'CCNC' :
                        if  NGS_info[j] in methy_info[7][i].split(','):
                            if methy_type == '6mA':
                                #methy_match_NGS[methy_index][j] = 'GM'
                                methy_match_NGS[methy_index][j] = 'AMMA'

                            else:
                                #methy_match_NGS[methy_index][j] = 'CM'
                                methy_match_NGS[methy_index][j] = 'CNNC'

                        else:
                            pass

        elif methy_loc == 'pos':
            ##pos
            methy_site = methy_info[1][i]
            if len(methy_match[(methy_match['Gene']==gene_name) & (methy_match['site']==methy_site)]) == 0:
                methy_match.loc[len(methy_match)] =  {'Gene':gene_name,'site':methy_site}
                tmp_methy = []
                for j in range(len(NGS_info)):
                    if NGS_info[j] in methy_info[7][i].split(','):
                        if methy_type == '6mA':
                            tmp_methy.append('AMAA')
                        else:
                            tmp_methy.append('CNCC')
                    else:
                        if methy_type == '6mA':
                            tmp_methy.append('AAAA')
                        else:
                            tmp_methy.append('CCCC')
                        
                methy_match_NGS.append(tmp_methy)
                
            else:
                methy_index = methy_match[(methy_match['Gene']==gene_name) & (methy_match['site']==methy_site)].index[0]
                for j in range(len(NGS_info)):
                    if methy_match_NGS[methy_index][j] == 'AAAA'  or methy_match_NGS[methy_index][j] == 'CCCC':
                        if  NGS_info[j] in methy_info[7][i].split(','):
                            if methy_type == '6mA':
                                methy_match_NGS[methy_index][j] = 'AMAA'
                            else:
                                methy_match_NGS[methy_index][j] = 'CNCC'
                        else:
                            pass
                    
        elif methy_loc == 'neg':
            ##pos
            methy_site = methy_info[1][i]
            if len(methy_match[(methy_match['Gene']==gene_name) & (methy_match['site']==methy_site)]) == 0:
                methy_match.loc[len(methy_match)] =  {'Gene':gene_name,'site':methy_site}
                tmp_methy = []
                for j in range(len(NGS_info)):
                    if NGS_info[j] in methy_info[7][i].split(','):
                        if methy_type == '6mA':
                            #tmp_methy.append('GM')
                            tmp_methy.append('AAMA')

                        else:
                            #tmp_methy.append('CM')
                            tmp_methy.append('CCNC')

                    else:
                        if methy_type == '6mA':
                            tmp_methy.append('AAAA')
                        else:
                            tmp_methy.append('CCCC')
                        
                methy_match_NGS.append(tmp_methy)
                
            else:
                methy_index = methy_match[(methy_match['Gene']==gene_name) & (methy_match['site']==methy_site)].index[0]
                for j in range(len(NGS_info)):
                    if methy_match_NGS[methy_index][j] == 'CCCC'  or methy_match_NGS[methy_index][j] == 'AAAA':
                        if  NGS_info[j] in methy_info[7][i].split(','):
                            if methy_type == '6mA':
                                # methy_match_NGS[methy_index][j] = 'GM'
                                methy_match_NGS[methy_index][j] = 'AAMA'

                            else:
                                # methy_match_NGS[methy_index][j] = 'CM'
                                methy_match_NGS[methy_index][j] = 'CCNC'

                        else:
                            pass

    return [methy_match,methy_match_NGS]


methy_match,methy_match_NGS = methy_info_match(methy_info,core_gene_info,NGS_info,methy_type=args.methy_name,methy_loc=args.methy_type,gene_type=args.gene_type)


methy_match_NGS = pd.DataFrame(methy_match_NGS)

###design map ped file
strain_fit = pd.read_csv(args.classA,header=None)


#ped
ped_info = []
for i in range(len(NGS_info)):
    tmp_ped = []
    if i+1 in strain_fit:
        fitness = 1
    else:
        fitness = 2
    
    tmp_ped.append(NGS_info[i])
    tmp_ped.append(NGS_info[i])
    tmp_ped.append(0)
    tmp_ped.append(0)
    tmp_ped.append(0)
    tmp_ped.append(fitness)
    
    ped_info.append(tmp_ped)

ped_out = args.out +'plink.ped'
with open(ped_out,'w') as f:
    for i in range(len(ped_info)):
        for j in range(len(ped_info[i])):
            f.write(str(ped_info[i][j])+' ')
        
        for j in range(len(methy_match_NGS)):
            f.write(str(methy_match_NGS[i][j][0])+' ')
            f.write(str(methy_match_NGS[i][j][1])+' ')
            f.write(str(methy_match_NGS[i][j][2])+' ')
            f.write(str(methy_match_NGS[i][j][3])+' ')
        f.write('\n')


map_info = []

for i in range(len(methy_match)):
    start_location = list(core_gene_info[core_gene_info['Gene'] == methy_match['Gene'][i]]['start'])[0]
    map_info.append([0,'snp'+str(i),0,round(start_location+methy_match['site'][i])])
    map_info.append([0,'snp'+str(i+len(methy_match)+1),0,round(start_location+methy_match['site'][i])+1])

map_info = pd.DataFrame(map_info)

map_out = args.out +'plink.map'
map_info.to_csv(map_out,header=None,index=None,sep='\t')























