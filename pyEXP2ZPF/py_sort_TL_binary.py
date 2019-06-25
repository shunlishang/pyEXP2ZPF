import collections
import os
import numpy as np
from bisect import bisect_left, insort
import json
import pandas as pd
import itertools
import copy

#########################################
stepT=50   # the temp step to select datapoint

#########################################
## define some funcitons below
##
def dedup(seq):
    # Remove duplicates. Preserve order first seen. Assume orderable, but not hashable elements'
    result = []
    seen = []
    for x in seq:
        i = bisect_left(seen, x)
        if i == len(seen) or seen[i] != x:
            seen.insert(i, x)
            result.append(x)
    return result

def extract_filename_metadata(fname):
    """
	File names have the format TL+237+BCC_A2+LIQUID+step2.json
	"""
    split_fname = fname.split('+')
    others     = split_fname[4]; #print('IN others=', others)
    phase2     = split_fname[3]; #print('IN 3=', phase2)
    phase1     = split_fname[2]; #print('IN 2=', phase1)
    block_no   = split_fname[1]; #print('IN 1=', block_no)
    property0  = split_fname[0]; #print('IN 0=', property0)
    phase_two  = [phase1, phase2]
    return phase_two, block_no

def To_write_json_file(in_dict, out_json_file):
    with open(out_json_file, 'w') as f:
        f.write('{\n')
        f.write('  "comment": ')
        f.write(json.dumps(in_dict.get('comment')))
        f.write(",")
        f.write('\n  "reference": ')
        f.write(json.dumps(in_dict.get('reference')))
        f.write(",")
        f.write('\n  "weight": ')
        f.write(json.dumps(in_dict.get('weight')))
        f.write(",")
        f.write('\n  "phases": ')
        f.write(json.dumps(in_dict.get('phases')))
        f.write(",")
        f.write('\n  "values": ')
        f.write(json.dumps(in_dict.get('values')))
        f.write(",")
        f.write('\n  "components": ')
        f.write(json.dumps(in_dict.get('components')))
        f.write(",")
        f.write('\n  "output": ')
        f.write(json.dumps(in_dict.get('output')))
        f.write(",")
        f.write('\n  "broadcast_conditions": false,')
        f.write('\n  "conditions": ')
        f.write(json.dumps(in_dict.get('conditions')))
        f.write('\n}')
        f.write('\n')
    f.close()

def To_write_exp_file(name_exp, temp2, val2a, val2b, temp3, val3a, val3b, val3c):
    with open(name_exp, 'w') as f:
        f.write('\n$DATAPLOT Phase diagram generated by py_sort_TL code')
        f.write('\nDATASET 1 ')
        f.write('\nATTRIBUTE CENTER ')
        f.write('\n')
        f.write('\nCOLOR 1 ')
        for i in range(len(temp2)):
            f.write('\n')
            asd1 = str(val2a[i]) + '  ' + str(temp2[i]) + '  S1'
            asd2 = str(val2b[i]) + '  ' + str(temp2[i]) + '  S1'
            f.write(asd1)
            f.write('\n')
            f.write(asd2)
        f.write('\n')
        f.write('\nCOLOR 2 ')
        for i in range(len(temp3)):
            f.write('\n')
            aa1 = str(val3a[i]) + '  ' + str(temp3[i]) + '  S100'
            aa2 = str(val3b[i]) + '  ' + str(temp3[i]) + '  S100'
            aa3 = str(val3c[i]) + '  ' + str(temp3[i]) + '  S100'
            f.write(aa1)
            f.write('\n')
            f.write(aa2)
            f.write('\n')
            f.write(aa3)
        f.write('\n')
    f.close()

### END of function
#############################################
#############################################
#############################################
print()
print('Need to run code pytcZPF...py first. This code is used to sort TL data via T step =', stepT, 'K')
current_directory_files = os.listdir('.')
two_phases_0 = []
block_0=[]
name0='nTL+'
name_exp='myExptTLINV.exp'
tot_temp2=[]; tot_val2a=[]; tot_val2b=[]

for fname in current_directory_files:
    if fname[0:3] == 'TL+':
        phase_two, block_no = extract_filename_metadata(fname)
        two_phases_0.append(sorted(phase_two))
        block_0.append(block_no)
        print('TL fname=', fname)
        #print('For two phases = ', phase_two)
#print('Total_phases=', len(two_phases_0), two_phases_0)
print('Total blocks=', block_0)
unique_two=dedup(two_phases_0)   # dedup is a function, see above
print('Unique list of TL phases= ', unique_two, len(unique_two))

for ii in unique_two:
    print()
    valu =[]; #nvalu=[]
    temp  =[]; #ntemp=[]
    print('For two TL phases = ', ii)
    for fname in current_directory_files:
        if fname[0:3] == 'TL+':
            phase_two, block_no = extract_filename_metadata(fname)
            if sorted(phase_two) == ii:
                with open(fname) as f:
                    data_dict = json.load(f)
                    valu.append(data_dict.get("values"))
                    temp.append(data_dict.get("conditions").get("T"))
                    text1=data_dict.get("comment")
                    text2=data_dict.get("reference")
                    phase_list=data_dict.get("phases")
                    all_elements=data_dict.get("components")
                    weightTL = data_dict.get("weight")
    name1=name0+ii[0]+'+'+ii[1]+'+stepT'+str(stepT)+'.json'
    ntemp = list(itertools.chain.from_iterable(temp))
    nvalu = list(itertools.chain.from_iterable(valu))
    aaa = copy.deepcopy(ntemp)
    print('Total TL aftre running the previous code pytcZPF...py =', len(aaa))
    #ttvv1 = pd.DataFrame({"a": ntemp, "b": nvalu})
    #ttvv2 = ttvv1.sort_values(by=['a'])
    #print('sorted ttvv by T=', ttvv2)
    ntemp=[ntemp[x] for x in np.argsort(aaa)]
    nvalu=[nvalu[x] for x in np.argsort(aaa)]
    mtemp =[ntemp[0]]
    mvalu =[nvalu[0]]
    #print()
    #print('ntemp_totalT=', ntemp)
    for i in range(1,len(ntemp)):
        if ntemp[i]-mtemp[-1] > stepT:
            mtemp.append(ntemp[i])
            mvalu.append(nvalu[i])
    print('Number of the selected temperatures = ', len(mtemp))
    #print(mvalu, len(mvalu))
    condits = dict(P=101325, T=mtemp)
    TL_dict = dict(comment=text1, phases=phase_list, reference=text2, values=mvalu, components=all_elements, \
                   output='ZPF', broadcast_conditions='false', conditions=condits, weight=weightTL)
    To_write_json_file(TL_dict, name1)
    ###------------
    tot_temp2.extend(mtemp)
    asd1=[]
    asd2=[]
    for i in range(len(mvalu)):
        asd1.append(mvalu[i][0][2][0])
        asd2.append(mvalu[i][1][2][0])
    tot_val2a.extend(asd1)
    tot_val2b.extend(asd2)
##
## for INV case below ###################################
with open('INV_Alls.json') as f:
    data_dict = json.load(f)
    temp3 = data_dict.get("conditions").get("T")
    val3a=[]; val3b=[]; val3c=[]
    for i in range(len(temp3)):
        val3a.append(data_dict.get("values")[i][0][2][0])
        val3b.append(data_dict.get("values")[i][1][2][0])
        val3c.append(data_dict.get("values")[i][2][2][0])

print()
print('temp3=',  len(temp3), '    ', temp3)
print('val3a=',  len(val3a), '    ', val3a)
print()
print('Len tot_temp2=',  len(tot_temp2))
print('Len tot_val2a=',  len(tot_val2a))
###########
To_write_exp_file(name_exp, tot_temp2, tot_val2a, tot_val2b, temp3, val3a, val3b, val3c)

###########

file0 = 'new_exp_data.json'
if os.path.exists(file0):
    with open(file0) as f:
        data_dict = json.load(f)
        temp4 = data_dict.get("conditions").get("T")
        vvv   = data_dict.get("values")
    val4a=[]; val4b=[]
    for i in range(len(temp4)):
        val4a.append(data_dict.get("values")[i][0][2][0])
        val4b.append(data_dict.get("values")[i][1][2][0])
    print()
    print('Expt temp4=', temp4)
    print('Expt val4a=', val4a)
    print('Expt val4b=', val4b)
    with open(name_exp, 'a+') as f:
        f.write('\n')
        f.write('\nCOLOR 3 ')
        for i in range(len(temp4)):
            if val4a[i] is not None:
                asd1 = str(val4a[i]) + '  ' + str(temp4[i]) + '  S3'
                f.write('\n')
                f.write(asd1)
            if val4b[i] is not None:
                asd2 = str(val4b[i]) + '  ' + str(temp4[i]) + '  S3'
                f.write('\n')
                f.write(asd2)
        f.write('\n')
    f.close()

###########
print()
print('========================================================')
print('Updated on 2019-05-29 to write EXP file')
print('Updated on 2019-05-20 for the bug of: ntemp[i]-mtemp[-1] ')
### THE END ##############

