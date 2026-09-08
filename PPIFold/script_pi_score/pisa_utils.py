# -*- coding: utf-8 -*-
#Adapted from pisa_utils.py (https://gitlab.com/topf-lab/pi_score/-/tree/master/score_scripts)

import json
import uuid
import subprocess

def run_pisa(pdbfile=None):
    assert pdbfile is not None
    
    session_id = uuid.uuid4().hex[:8]
    outfle = pdbfile.rsplit('.pdb',1)[0] + '_pisa_' + session_id + '.xml'
    
    # lancer pisa avec session unique
    cmds = [
        'pisa %s -analyse %s' % (session_id, pdbfile),
        'pisa %s -xml interfaces > %s' % (session_id, outfle),
        'pisa %s -erase' % (session_id)
    ]

    for i, cmd in enumerate(cmds):
        proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT,close_fds=True)
        proc.wait() 
    return outfle


def parse_pisa_xml_outfle(xml_fle,outfle_name=None):
    import xml.etree.ElementTree as ET
    intf_area = ''
    pvalue = ''
    energy = ''
    hb = ''
    sb = ''
    dict_pisa = {}
    if outfle_name == None:
        outfle_name = xml_fle.split('_pisa')[0] + '_dict_pisa.json'
    try:
        root = ET.parse(xml_fle).getroot()
        for interface in root.findall('interface'):
            intf_ch = []
            intf_area = "%.2f" % float(interface.find('int_area').text)
            pvalue = "%.2f" % float(interface.find('pvalue').text)
            energy = "%.2f" % float(interface.find('int_solv_en').text)
            for mol in interface.findall('molecule'):
                mol1 = mol.find('chain_id').text
                intf_ch.append(mol1)
            intf_chains = '_'.join(intf_ch)
            dict_pisa[intf_chains] = {}
            for hbonds in interface.findall('h-bonds'):
                hb = int(hbonds.find('n_bonds').text)
            for sb in interface.findall('salt-bridges'):
                sb = int(sb.find('n_bonds').text)
            dict_pisa[intf_chains]['int_area'] = intf_area
            dict_pisa[intf_chains]['pvalue'] = pvalue
            dict_pisa[intf_chains]['int_solv_en'] = energy
            dict_pisa[intf_chains]['hb'] = hb
            dict_pisa[intf_chains]['sb'] = sb
    except:
        root = []
        pass
    if dict_pisa:
        with open(outfle_name,'w') as outfle:
            json.dump(dict_pisa,outfle)
    return dict_pisa
