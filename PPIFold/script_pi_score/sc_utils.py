# -*- coding: utf-8 -*-
#Adapted from sc_utils.py (https://gitlab.com/topf-lab/pi_score/-/tree/master/score_scripts)
import os
import json
import uuid
import subprocess

def write_sc_script(pdbfile,
                    chain_lst=None,
                    outfle=None,
                    ):
    '''
    chain_lst should be list of only two chain Ids e.g. [A,B]
    '''
    #print ('Writing sc script', chain_lst)
    assert (chain_lst != None),"Input chain ID list!"
    assert len(chain_lst) == 2, "Input list with only two chain IDs e.g [A,B]"
    if outfle == None:
        outfle = pdbfile + '_' + '_'.join(chain_lst) + '_sc_infle.sh'
    #print (outfle, 'OUTFLE for SC')
    out = open(outfle,'w')
    out.write('#!/bin/bash\n')
    out.write('sc XYZIN ' + pdbfile + ' <<eof \n')
    out.write('MOLECULE 1\nCHAIN ' + chain_lst[0] + '\n' + 'MOLECULE 2\nCHAIN ' + chain_lst[1] + '\nEND\neof')
    out.close()
    os.system('chmod 777 ' + outfle)
    return outfle

def run_sc(scriptfile):
    session_id = uuid.uuid4().hex[:8]
    unique_script = "%s_%s.sh" % (scriptfile.rsplit('.', 1)[0], session_id)
    tmp_out = "%s_tmp_sc_%s.out" % (scriptfile.rsplit('.', 1)[0], session_id)

    subprocess.call(["cp", scriptfile, unique_script])
    subprocess.call(["sed", "-i", "s/session name/session_%s/g" % session_id, unique_script])



    with open(tmp_out, "w") as out_fh:
        proc = subprocess.Popen([unique_script], stdout=out_fh, stderr=subprocess.STDOUT, close_fds=True)
        proc.wait()

    os.remove(unique_script)
    parse_sc_output(tmp_out)

def parse_sc_output(outfle):
    new_name = outfle.rsplit('_tmp_sc.out',1)[0] + '_dict_sc_scores.json'
    if os.path.isfile(new_name):
        with open(new_name) as json_file:
            dict_out = json.load(json_file)
    #else:
    dict_out = {}
    ch_lst = []
    pdb = ''
    sc = ''
    #print(outfle)
    with open(outfle,'r') as out:
        for lne in out:
            if 'Number of atoms selected in chain' in lne:
                ch_lst.append(lne.split('=')[0].split()[-1])
            if 'Logical name: XYZIN  File name: ' in lne:
                #pdb = lne.split()[-1][0:4]
                #for em challenge targets
                pdb = lne.split()[-1].split('_clean')[0]

            if 'Shape complementarity statistic Sc' in lne:
                sc = lne.split('=')[-1].strip()
        try:
            if pdb in dict_out.keys():
                dict_out[pdb].update({'_'.join(ch_lst):sc})
            else:
                dict_out[pdb] = {'_'.join(ch_lst):sc}
        except:
            pass

        with open(new_name,'w') as sec_outfle:
            json.dump(dict_out,sec_outfle)




