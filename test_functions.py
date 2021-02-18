"""
Tests for various routines
"""
from processing import sub_replacement, change_support_format, trim_name
import re


def test_node_processing():
    node_re = re.compile(':([\d.E-]+)\[([\d.E-]+)\]')
    match = node_re.match(':0.95911[100]')
    assert sub_replacement(match) == '100:0.95911'
    e_match = node_re.match(':2.5E-2[3.3E-10]')
    assert sub_replacement(e_match) == '3.3E-10:2.5E-2'


def test_support_reformating():
    simple_tree = '((A:1, B:0.7):0.8[65], C)'
    assert change_support_format(simple_tree) == '((A:1, B:0.7)65:0.8, C)'
    escaped_tree = '((\'A\':1, \'B\':0.7):0.8[65], \'C\')'
    assert change_support_format(escaped_tree) == '((\'A\':1, \'B\':0.7)65:0.8, \'C\')'
    decimal_tree = '((A:1, B:0.7):0.8[0.7331], C)'
    assert change_support_format(decimal_tree) == '((A:1, B:0.7)0.7331:0.8, C)'


def test_name_trimming():
    # Domain position is trimmed
    name = 'Nitzschia_punctata,_Strain_CCMP561|CAMPEP_0199315576_(5-177)_1'
    assert trim_name(name) == 'Nitzschia_punctata,_Strain_CCMP561|CAMPEP_0199315576_1'
    # Clean sequences are unaffected
    clean_name = 'Skeletonema_costatum,_Strain_1716|CAMPEP_0113383910_2'
    assert trim_name(clean_name) == clean_name
