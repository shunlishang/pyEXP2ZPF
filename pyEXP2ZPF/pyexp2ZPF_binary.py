import re
import numpy as np
import json
###################################################################################
first_elem ='AL'          # can be upper or lower letters for elements, will adjust them to upper later
second_elem='W'          # This element is x(elem) to plot phase diagram
weightTL   = 5            # Weight of tie lines
weightINV  = 10           # Weight of Inv equilibria
nTL_step   = 1            # for tie lines, the step: 1 for all; 2 half TL, 3 for 1/3 TL, ...

my_file=open("tcexp.exp",'r')     # open Thermo-Calc EXP file: binary system

text1 = 'Comment: from TDB file'
text2 = 'Reference: from TDB'

#---------------------------------------
first_elem =first_elem.upper()
second_elem=second_elem.upper()
all_elements=np.sort([first_elem, second_elem])            # alpha-beta list of all upper elements
all_elements=[x for x in all_elements]
textp  = my_file.read()                  # read Thermo-Calc EXP file
####################################################################################
####################################################################################
################ Begin: def some funcitons #########################################
def find_TL_in_PLOT(blockTL_plot):
    dd=dict()
    Comp1=[]; nComp1=[]; Temp2=[]; nTemp2=[]; phase=[]
    index_start = blockTL_plot.find("   M")              # seperate by "   M"
    tieline=blockTL_plot[index_start - 36:]
    tt_lines = tieline.split('\n');  # print('tt_lines =', tt_lines, "and lines = ", len(tt_lines))
    for line in tt_lines:
        if line == '': continue
        filtered_line = [val for val in line.split(' ') if val != '']
        Comp1.append(filtered_line[0])          # yields a list of composition
        nComp1.append(float(filtered_line[0]))  # yields a list of composition, number
        Temp2.append(filtered_line[1])          # yields a list of composition
        nTemp2.append(float(filtered_line[1]))  # yields a list of composition, number
    phase_start = blockTL_plot.find("(")
    phase_end = blockTL_plot.find(")")
    phelem = blockTL_plot[phase_start + 1:phase_end];  # print(phelem)
    asdf = re.split(",", phelem)
    phase = asdf[0];  # print(phase)
    if '#' in phase:
        aaa = re.split("#", phase)
        phase=aaa[0]
    dd['TLphase']=phase
    dd['TLcomp1']=nComp1
    dd['TLtemp2']=nTemp2
    return dd
#---------------

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
                phase = bb[0]
                phase_find.append(phase)
            else:
                phase_find.append(aa)
    return phase_find
##----------------

def two_values_in_one_line(line_value):
    filtered_line = [val for val in line_value.split(' ') if val != '']
    value1=float(filtered_line[0])  # yields the first value
    value2=float(filtered_line[1])  # yields the second value
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

#---------------
def To_build_TLdict(ALLTL_PH, ALLTL_comp1, ALLTL_temp2, all_elements, second_elem, nTL_step, text1, text2, weightTL):
    # for tie-line case in a small block, NOT for all tie-lines
    phase_list = np.unique(ALLTL_PH)
    phase_list = [x for x in phase_list]
    val_temp = [];  val_TL = []
    for nnn in range(0, len(ALLTL_comp1[0]), nTL_step):
        qqq = [[ALLTL_PH[0], [second_elem], [ALLTL_comp1[0][nnn]]],
               [ALLTL_PH[1], [second_elem], [ALLTL_comp1[1][nnn]]]]
        val_TL.append(qqq)
        val_temp.append([ALLTL_temp2[0][nnn]])
    val_temp = [x[0] for x in val_temp]
    condits = dict(P=101325, T=val_temp)
    TL_dict = dict(comment=text1, phases=phase_list, reference=text2, values=val_TL, components=all_elements,
                   output='ZPF', broadcast_conditions='false', conditions=condits, weight=weightTL)
    return TL_dict

################## End: def some funcitons #########################################
####################################################################################
####################################################################################
#
# next to split the blocks and process both Tie-Line and Inv Eqs separately
#
blocksp = re.split("BLOCKEND", textp)        # splits txt file by BLOCKEND
blocksp = blocksp[0:len(blocksp)-1]          # remove the last block
index_TL=[]; index_INV=[]
for index in range(0,len(blocksp)):  #for each block
    if "INVARIANT EQUILIBRIUM" in blocksp[index]:  ## case 1 of 2 for Inv blocks
        #print("The block of", index, " is for invariant Eqs")
        num_M=blocksp[index].count("   M")
        if num_M == 3: index_INV.append(index)
    else: #-------- case 2 of 2: for Tie Line block/case below
        #print("The block of", index, " is for tie lines")
        index_TL.append(index)
print('The blocks for index_INV = ', index_INV)
print('The blocks for index_TL  = ', index_TL)
print()
#################################################################################
ALLINV_comp1 =[]; ALLINV_temp2 =[];  ALLINV_PH=[]
for index in index_INV:  #for each block of INV
    print("The block of", index, " is for invariant Eqs")
    blockINV = re.split("\$ INVARIANT EQUILIBRIUM", blocksp[index])

    ph3_comp1 = [];  ph3_temp2 = [];  ph3_phases = [];   nph3_comp1 = []; nph3_temp2 = []
    aa_lines = blockINV[0].split('\n');  #print('blockINV0=', blockINV[0], 'aa_lines Inv = ', len(aa_lines))
    for line in aa_lines:
        if " BLOCK " in line: continue
        if line == '': continue
        filtered_line = [val for val in line.split(' ') if val != '']
        phph=filtered_line[-1]
        if '#' in phph:
            bb = re.split("#", phph)
            ph3_phases.append(bb[0])
        else:
            ph3_phases.append(phph)
    #print('3-phases_inv_block = ', ph3_phases)

    aa_lines = blockINV[1].split('\n');  #print('blockINV1=', blockINV[1], 'aa_lines Inv = ', len(aa_lines))
    list_del = ("BLOCK ", "COLOR ", "   M")
    my_phase=[]
    for line in aa_lines:
        if any(s in line for s in list_del): continue
        if line == '': continue
        check_line = line; #print('Here the line value =', line)
        asdf=two_values_in_one_line(check_line)
        ph3_comp1.append(asdf[0])
        ph3_temp2.append(asdf[1])
        for jj in index_TL:
            aa=find_phase_in_plot(check_line, blocksp[jj]); #print('The phase find is = ', aa)
            if len(aa) == 1: my_phase.append(aa); break
    print('Here are the 3 Phases for INV =', my_phase)
    print('Here are the 3 Comp1  for INV =', ph3_comp1)
    print('Here are the 3 Temp2  for INV =', ph3_temp2)
    ALLINV_comp1.append(ph3_comp1)
    ALLINV_temp2.append(ph3_temp2)
    my_phase=[x[0] for x in my_phase]
    ALLINV_PH.append(my_phase)
    print()
#print('ALL Inv comp1 =', ALLINV_comp1); print('ALL Inv comp2 =', ALLINV_temp2); print('ALL Inv Phases=', ALLINV_PH)
#--------------------------
phase_list3=np.unique(ALLINV_PH)
phase_list3=[x for x in phase_list3]
print('List All Unique Phases in Inv =', phase_list3)

val_temp3 = []; val_INV = []
for iINV in range(0,len(ALLINV_PH)):
    qqq = [[ALLINV_PH[iINV][0], [second_elem], [ALLINV_comp1[iINV][0]]],
           [ALLINV_PH[iINV][1], [second_elem], [ALLINV_comp1[iINV][1]]],
           [ALLINV_PH[iINV][2], [second_elem], [ALLINV_comp1[iINV][2]]]]
    val_INV.append(qqq)
    val_temp3.append(ALLINV_temp2[iINV][0])
condit3=dict(P=101325, T=val_temp3)
#print('T and P = ', condit3); print('val_TL3= ', val_TL3)
INV_dict=dict(comment=text1, phases=phase_list3, reference=text2, values=val_INV, components=all_elements,
              output='ZPF', broadcast_conditions='false', conditions=condit3, weight=weightINV)
To_write_json_file(INV_dict, 'INV_Alls.json')

print('########################## END of the INV cases #########################')
print()
####################################################################################
####################################################################################
###########################################: next for TL: one by one and all_together
ALLTL_comp1 = [];     ALLTL_temp2 = [];     ALLTL_PH = []

for index in index_TL: #[0, 2]:  #for each block of Tie-Lines
    print("The block of", index, " is for tie lines")
    blockTL = re.split("\$ PLOTTED COLUMNS ARE :", blocksp[index])
    index_0    = blockTL[0].find(" $ BLOCK #"); #print('index_0=',index_0)
    newcontent = blockTL[0][index_0 + 10 : ]; #print('newcontent=',newcontent)
    block_number = re.split(" ", newcontent)[0]
    print('block_number=',block_number)

    phaseTL=[]; compTL_1=[]; tempTL_2=[]
    for ii in [1,2]:
        #print('for case in block', ii) # print('blockTL= \n', blockTL[ii])
        TLdata=find_TL_in_PLOT(blockTL[ii])
        #print('blockTL of phase= ',TLdata.get('TLphase'))
        phaseTL.append(TLdata.get('TLphase'))
        compTL_1.append(TLdata.get('TLcomp1'))
        tempTL_2.append(TLdata.get('TLtemp2'))
    print('These two phases in TL =', phaseTL);
    #print('comp1 in TL for two phases len= ',len(compTL_1), ' values=', compTL_1); print('temp2 in TL for two phases=', tempTL_2)

    TLname=['TL+',block_number,'+',phaseTL[0],'+',phaseTL[1],'+step', str(nTL_step),'.json']
    TLname=''.join(TLname); #print('TL_NAME=', TLname)
    TL_dict=To_build_TLdict(phaseTL, compTL_1, tempTL_2, all_elements, second_elem, nTL_step, text1, text2, weightTL)
    To_write_json_file(TL_dict, TLname)
    #-----------
    ALLTL_comp1.append(compTL_1)
    ALLTL_temp2.append(tempTL_2)
    ALLTL_PH.append(phaseTL)
    print()
#-----------------------------------------------------------------
#print('ALLTL_PH len and values=', len(ALLTL_PH), ALLTL_PH)
phase_list=np.unique(ALLTL_PH)
phase_list=[x for x in phase_list]
print('List ALL Unique Phases of TLs: ', phase_list)
#print('ALLTL_comp1 len and values=', len(ALLTL_comp1), ALLTL_comp1)

val_temp = []; val_TL = []
for iTL in range(0,len(ALLTL_PH)):
    for nnn in range(0, len(ALLTL_comp1[iTL][0]), nTL_step):
        qqq=[[ALLTL_PH[iTL][0], [second_elem], [ALLTL_comp1[iTL][0][nnn]]],
             [ALLTL_PH[iTL][1], [second_elem], [ALLTL_comp1[iTL][1][nnn]]]]
        val_TL.append(qqq)
        val_temp.append(ALLTL_temp2[iTL][0][nnn])
condit=dict(P=101325, T=val_temp)
INV_TL=dict(comment=text1, phases=phase_list3, reference=text2, values=val_TL, components=all_elements,
              output='ZPF', broadcast_conditions='false', conditions=condit, weight=weightTL)
name2=['TL_ALLs_step',str(nTL_step),'.json']
ALLname=''.join(name2)
To_write_json_file(INV_TL, ALLname)

