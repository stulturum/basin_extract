'''
Created on 31 Oct 2020

@author: thomasgumbricht
'''

import xml.etree.ElementTree as ET
from collections import defaultdict

def Etree_to_dict(t):
    ''' Read XML file using xml.etree.ElementTree
    '''
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(Etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v
                     for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v)
                        for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d

def ReadXML(xmlFPN):
    '''  Call XML Reader and return structured params
    '''
    tree = ET.parse(xmlFPN)
    
    d = Etree_to_dict( tree.getroot() )
    
    return d