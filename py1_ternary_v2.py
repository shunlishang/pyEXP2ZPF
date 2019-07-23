import re
import glob
import numpy as np
import json
#########################################################################################
first_element='AL'        # element for the first corner
weightTL   = 1            # Weight of tie line
weightINV  = 1            # Weight of Inv equilibria

my_va = 'VA'              # 'VA' if VA included, otherwise '' only

nTL_step   = 2            # 1 for all TL; 2 half TL, 3 for 1/3 TL, ...

iadjust  = -1   #   > 0 to adjust the value; < 0 NOT

temp_fixed_input = 1373    # the temperature to calculate isothermal section

text1 = 'Comt NIST AlCoW TDB file'
text2 = 'Ref Peisheng Wang CALPHAD paper'

################ Begin: def some functions ########################################
################ Begin: def some functions ########################################
###################################################################################
##
def find_TL_in_PLOT(blockTL_plot, iadjust):
    dd=dict()
    Comp1=[]; nComp1=[]; Comp2=[]; nComp2=[]; phase=[]
    index_start = blockTL_plot.find("   M")             # seperate by "   M"
    if index_start > 1:
        #print('##########index start =', index_start)
        tieline = blockTL_plot[index_start - 36:]
        tt_lines = tieline.split('\n');  # print('tt_lines =', tt_lines, "and lines = ", len(tt_lines))
        for line in tt_lines:
            if line == '': continue
            filtered_line = [val for val in line.split(' ') if val != '']
            Comp1.append(filtered_line[0])  # yields a list of composition
            #print('##########in find_TL_in_PLOT=', filtered_line[0])
            value1 = float(filtered_line[0])
            if value1 < 1e-5 and iadjust > 0:
                print('find TL value < 1e-5 with its value = ', value1)
                value1 = 1e-5
            if value1 > 0.99999 and iadjust > 0:
                print('find TL value > 0.99999 with its value = ', value1)
                value1 = 0.99999
            nComp1.append(value1)  # yields a list of composition, number
            ##
            Comp2.append(filtered_line[1])  # yields a list of composition
            value2 = float(filtered_line[1])
            if value2 < 1e-5 and iadjust > 0:
                print('find TL value < 1e-5 with its value = ', value2)
                value2 = 1e-5
            if value2 > 0.99999 and iadjust > 0:
                print('find TL value > 0.99999 with its value = ', value2)
                value2 = 0.99999
            if (value1 + value2) > 1:
                value2 = 1 - value1
                print('****Place 1: value1 + value2 > 1: ', value1, value2)
            nComp2.append(value2)  # yields a list of composition, number
        phase_start = blockTL_plot.find("(")
        phase_end = blockTL_plot.find(")")
        phelem = blockTL_plot[phase_start + 1:phase_end];  # print(phelem)
        asdf = re.split(",", phelem)
        phase = asdf[0];  # print(phase)
        if '#' in phase:
            bb = re.split("#", phase)
            phase = bb[0]
    dd['TLphase']=phase
    dd['TLcomp1']=nComp1
    dd['TLcomp2']=nComp2
    return dd
#---------------

def two_values_in_one_line(line_value, iadjust):
    filtered_line = [val for val in line_value.split(' ') if val != '']
    value1=float(filtered_line[0])  # yields the first value
    if value1 < 1e-5 and iadjust > 0:
        print('find value < 1e-5 with its value = ', value1)
        value1 = 1e-5
    if value1 > 0.99999 and iadjust > 0:
        print('find value > 0.99999 with its value = ', value1)
        value1 = 0.99999

    value2=float(filtered_line[1])  # yields the second value
    if value2 < 1e-5 and iadjust > 0:
        print('find value < 1e-5 with its value = ', value2)
        value2 = 1e-5
    if value2 > 0.99999 and iadjust > 0:
        print('find value > 0.99999 with its value = ', value2)
        value2 = 0.99999

    if (value1 + value2) > 1:
        value2 = 1 - value1
        print('****Place 2: value1 + value2 > 1: ', value1, value2)
    dd=[value1, value2]
    return dd
##----------------

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
##----------------

def find_phase_in_plot(line_value, whole_blockTL):
    blockTL = re.split("\$ PLOTTED COLUMNS ARE", whole_blockTL)
    phase_find=[]
    for ii in range(1, len(blockTL)):
        if line_value in blockTL[ii]:
            #print('In func, for block=',ii, 'values are=', blockTL[ii])
            asdf=blockTL[ii]
            phase_start = asdf.find("(")
            phase_end = asdf.find(")")
            phelem = asdf[phase_start + 1:phase_end]
            aa = re.split(",", phelem)[0]
            if '#' in aa:
                bb = re.split("#", aa)
                phase_find.append(bb[0])
            else:
                phase_find.append(aa)
    return phase_find
##----------------

def find_two_elems(textp):
    aa_lines = textp.split('\n')
    elem_tot=[]
    for line in aa_lines:
        if " XTEXT " in line:
            filtered_line = [val for val in line.split(' ') if val != '']
            elem1=filtered_line[2]
        if " YTEXT " in line:
            filtered_line = [val for val in line.split(' ') if val != '']
            elem2=filtered_line[2]
    elem_tot=[elem1, elem2]
    return elem_tot

################# End: def some funcitons #########################################
###################################################################################
###################################################################################
current_directory_exp = glob.glob('*.exp')
exp_tc_name =current_directory_exp[0]
print('Used EXP file for this job is: ', current_directory_exp[0])
my_file=open(exp_tc_name,'r')    # open Thermo-Calc EXP file: binary system
### ---------------------

with open('temp1.txt', 'w') as f:
    f.write(str(3))  # 3 for ternary and 2 for binary
    f.close()
##---------------------------------------------------------------------------------------
textp  = my_file.read()                      # read Thermo-Calc EXP file
elemtwo=find_two_elems(textp)
print('Two rest elements =', elemtwo)
elem1_change =elemtwo[0]
elem2_change =elemtwo[1]
nfirst_elem=first_element.upper()
nelem1_change=elem1_change.upper()
nelem2_change=elem2_change.upper()
if my_va:
    all_elements = np.sort([nfirst_elem, nelem1_change, nelem2_change, my_va])
else:
    all_elements = np.sort([nfirst_elem, nelem1_change, nelem2_change])
all_elements=[x for x in all_elements]

blocksp = re.split("BLOCKEND", textp)        # splits txt file by BLOCKEND
blocksp = blocksp[0:len(blocksp)-1]          # remove the last block

index_TL=[]; index_INV=[]
for index in range(0,len(blocksp)):  #for each block
    if "INVARIANT EQUILIBRIUM" in blocksp[index]: # case 1 of 2: for the Inv Eq block
        #print("The block of", index, " is for invariant Eqs")
        num_M=blocksp[index].count("   M")
        if num_M == 3: index_INV.append(index)
    else: #-------- case 2 of 2: for Tie Line block/case below
        #print("The block of", index, " is for tie lines")
        index_TL.append(index)
print('The blocks for index_INV = ', index_INV)
print('The blocks for index_TL  = ', index_TL)
print()

#########################################################################################
### next for INV case

ALLINV_comp1 =[]; ALLINV_comp2 =[];  ALLINV_PH=[]
for index in index_INV:  #for each block of INV
    blockINV = re.split("\$ INVARIANT EQUILIBRIUM", blocksp[index])
    index_0    = blockINV[0].find(" $ BLOCK #"); #print('index_0=',index_0)
    newcontent = blockINV[0][index_0 + 10 : ]; #print('newcontent=',newcontent)
    block_number = re.split(" ", newcontent)[0]
    print("The block of", block_number, " is for INV")

    ph3_comp1 = [];   ph3_comp2 = [];  ph3_phases = []
    nph3_comp1 = []; nph3_comp2 = []

    aa_lines = blockINV[0].split('\n');  #print('blockINV0=', blockINV[0], 'aa_lines Inv = ', len(aa_lines))
    for line in aa_lines:
        if " BLOCK " in line: continue
        if line == '': continue
        filtered_line = [val for val in line.split(' ') if val != '']
        ph3_phases.append(filtered_line[-1])
    print('3-phases in inv_block = ', ph3_phases)

    aa_lines = blockINV[1].split('\n');  #print('blockINV1=', blockINV[1], 'aa_lines Inv = ', len(aa_lines))
    list_del = ("BLOCK ", "COLOR ", "   M")
    my_phase=[]
    for line in aa_lines:
        if any(s in line for s in list_del): continue
        if line == '': continue
        check_line = line; #print('Here the line value =', line)
        asdf=two_values_in_one_line(check_line, iadjust)
        ph3_comp1.append(asdf[0])
        ph3_comp2.append(asdf[1])

        for jj in index_TL:
            aa=find_phase_in_plot(check_line, blocksp[jj]); #print('The phase find is = ', aa)
            if len(aa) == 1: my_phase.append(aa); break
    print('Here is the phase find for INV =', my_phase)
    print('Here is the comp1 find for INV =', ph3_comp1)
    print('Here is the comp2 find for INV =', ph3_comp2)
    print()
    ALLINV_comp1.append(ph3_comp1)
    ALLINV_comp2.append(ph3_comp2)
    my_phase = [x[0] for x in my_phase]
    ALLINV_PH.append(my_phase)
    phase_list30 = np.unique(my_phase)
    phase_list30 = [x for x in phase_list30]
    print('List current unique Phases Inv: ', phase_list30)
    val_INV0 = [[my_phase[0], [nelem1_change, nelem2_change], [ph3_comp1[0], ph3_comp2[0]]],
                [my_phase[1], [nelem1_change, nelem2_change], [ph3_comp1[1], ph3_comp2[1]]],
                [my_phase[2], [nelem1_change, nelem2_change], [ph3_comp1[2], ph3_comp2[2]]]]
    condit30  = dict(P=101325, T=temp_fixed_input)
    INV_dict0 = dict(comment=text1, phases=phase_list30, reference=text2, values=[val_INV0], components=all_elements, \
                    output='ZPF', broadcast_conditions='false', conditions=condit30, weight=weightINV)
    name_inv_here = 'INV+' + str(block_number) + '+' + my_phase[0] + '+' + my_phase[1] + \
            '+' + my_phase[2] + '.json'
    To_write_json_file(INV_dict0, name_inv_here)
#print('ALL Inv comp1 =', ALLINV_comp1); print('ALL Inv comp2 =', ALLINV_comp2); print('ALL Inv Phases=', ALLINV_PH)
#-------------------
phase_list3=np.unique(ALLINV_PH)
phase_list3=[x for x in phase_list3]
print('List ALL Phases Inv: ', phase_list3)

val_temp3 = []; val_INV = []
for iINV in range(0,len(ALLINV_PH)):
    qqq = [[ALLINV_PH[iINV][0], [nelem1_change, nelem2_change], [ALLINV_comp1[iINV][0], ALLINV_comp2[iINV][0]]],
           [ALLINV_PH[iINV][1], [nelem1_change, nelem2_change], [ALLINV_comp1[iINV][1], ALLINV_comp2[iINV][1]]],
           [ALLINV_PH[iINV][2], [nelem1_change, nelem2_change], [ALLINV_comp1[iINV][2], ALLINV_comp2[iINV][2]]]]
    val_INV.append(qqq)
    val_temp3.append(temp_fixed_input)
condit3=dict(P=101325, T=val_temp3)
#print('T and P = ', condit3); print('val_TL3= ', val_TL3)

INV_dict=dict(comment=text1, phases=phase_list3, reference=text2, values=val_INV, components=all_elements,
              output='ZPF', broadcast_conditions='false', conditions=condit3, weight=weightINV)
To_write_json_file(INV_dict, 'INV_Alls.json')

#####################################################################################################
#####################################################################################################
############################################################### next for TL
ALLTL_comp1 = [];     ALLTL_comp2 = [];     ALLTL_PH = []

for index in index_TL: #[0, 2]:  #for each block of Tie-Lines
    print()
    blockTL = re.split("\$ PLOTTED COLUMNS ARE :", blocksp[index])
    index_0    = blockTL[0].find(" $ BLOCK #"); #print('index_0=',index_0)
    newcontent = blockTL[0][index_0 + 10 : ]; #print('newcontent=',newcontent)
    block_number = re.split(" ", newcontent)[0]
    print("The block of", block_number, " is for tie lines")
    testdata = find_TL_in_PLOT(blockTL[1], iadjust)
    if not testdata.get('TLphase'):
        print("#### ========== NOTE: block of", block_number, " is EMPTY")
        pass
    else:
        phaseTL = [];  compTL_1 = [];  compTL_2 = []
        for ii in [1, 2]:
            # print('for case in block', ii) # print('blockTL= \n', blockTL[ii])
            TLdata = find_TL_in_PLOT(blockTL[ii], iadjust); #print('blockTL of phase= ', TLdata.get('TLphase'))
            phaseTL.append(TLdata.get('TLphase'))
            compTL_1.append(TLdata.get('TLcomp1'))
            compTL_2.append(TLdata.get('TLcomp2'))
        print('***** Two phases in TL =', phaseTL);
        # print('comp1 in TL for two phases=', compTL_1) print('comp2 in TL for two phases=', compTL_2)
        ALLTL_comp1.append(compTL_1)
        ALLTL_comp2.append(compTL_2)
        ALLTL_PH.append(phaseTL)
        temp1 = [];  val_temp = [];  val_TL = [];  numberT = []
        for nnn in range(0, len(compTL_1[0]), nTL_step):
            qqq = [[phaseTL[0], [nelem1_change, nelem2_change], [compTL_1[0][nnn], compTL_2[0][nnn]]],
                   [phaseTL[1], [nelem1_change, nelem2_change], [compTL_1[1][nnn], compTL_2[1][nnn]]]]
            val_TL.append(qqq)
            temp1.append(temp_fixed_input)
            val_temp.append(temp_fixed_input)
        condit = dict(P=101325, T=val_temp)
        phase22 = np.unique(phaseTL);
        phase22 = [x for x in phase22]
        TL_dict = dict(comment=text1, phases=phase22, reference=text2, values=val_TL, components=all_elements, \
                       output='ZPF', broadcast_conditions='false', conditions=condit, weight=weightTL)
        TLname = TLname = 'TL+' + str(block_number) + '+' + phaseTL[0] + '+' + phaseTL[1] + '+step' + \
                          str(nTL_step) + '.json'
        To_write_json_file(TL_dict, TLname)

#-------------------
phase_list=np.unique(ALLTL_PH)
phase_list=[x for x in phase_list]
print('List unique phases TL case = ', phase_list)
print()
print('--------------------------------------------------------------------')
print('2019-07-22: judge empty TL block and write each INV for a file')
print('2019-07-15: remove #2 etc in FCC_A1#2 etc, read 2 elements in EXP file, adjust compositions')

