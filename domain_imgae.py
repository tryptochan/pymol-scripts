#!/usr/bin/env python

"""Generate imgaes for domains and PDB file for domain in assembly.

Domain generation script used in ECOD.

Required inputs: 9 digit domain UID, PDB ID, PDB range.
Optional inputs: ligands as "Chain ID:Residue number,...", like A:801,A:802
                 Other PDB ranges of domains in the same domain assembly

Requires Pymol 1.7.4 and later to handle multi-letter chain ID correctly.
"""

import sys
import os
import gzip
import re


args = sys.argv[1:]
if len(args) < 3:
    print 'Usage: pymol -qc domain_image.py -- uid pdb range [ligands, assembly partner ranges]'
    sys.exit()
from rainbow import rainbow
import pymol
from pymol import cmd


# Use PDB first, then try CIF. Could only use CIF.
# Now only use mmCIF for large structures.
pdbpath = '/usr2/pdb/data/structures/divided/pdb/%s/pdb%s.ent.gz'
cifpath = '/usr2/pdb/data/structures/divided/mmCIF/%s/%s.cif.gz'

basepath = ''

ptn = re.compile(r'(-?[0-9]+[A-Z]?)-(-?[0-9]+[A-Z]?)')

#HET ligands with only one atom (Some are obsolete).
#list generated from /pdb/data/monomers/het_dictionary.txt
#Hopefully these ligands will not need to be added
atomic_ligands = set(['3CO', '3NI', '4MO', '6MO', 'AG', 'AL', 'AR', 'ARS',
        'AU', 'AU3', 'BA', 'BR', 'BRO', 'BS3', 'CA', 'CD', 'CE', 'CL',
        'CLO', 'CO', 'CR', 'CS', 'CU', 'CU1', 'CU3', 'D8U', 'DUM', 'DY',
        'ER3', 'EU', 'EU3', 'F', 'FE', 'FE2', 'FLO', 'GA', 'GD', 'GD3',
        'H', 'HG', 'HO', 'HO3', 'IDO', 'IN', 'IOD', 'IR', 'IR3', 'K',
        'KR', 'LA', 'LI', 'LU', 'MG', 'MN', 'MN3', 'MO', 'NA', 'NGN',
        'NI', 'O', 'OS', 'OS4', 'OX', 'OXO', 'PB', 'PD', 'PR', 'PT',
        'PT4', 'QTR', 'RB', 'RE', 'RH', 'RH3', 'RU', 'S', 'SB', 'SE',
        'SM', 'SR', 'TA0', 'TB', 'TE', 'TL', 'U1', 'UNX', 'V', 'W',
        'XE', 'Y1', 'YB', 'YB2', 'YT3', 'ZN', 'ZN2'])



def range2sele(range_str):
    """Convert range string to Pymol selection."""
    sele = []
    chains = []
    for i in range_str.split(','):
        c, r = i.split(':')
        # minus index need to be excaped by \
        m = ptn.match(r)
        if m is None:
            print 'WARNING: range not correctly formatted, %s' % i
        else:
            start, end = m.groups()
            if start.startswith('-'):
                start = '\\' + start
            if end.startswith('-'):
                end = '\\' + end
            r = "%s-%s" % (start, end)
        sele.append('c. %s & i. %s' % (c, r))
        chain = c.strip()
        if chain not in chains:
            chains.append(chain)
    return (' | '.join(sele), chains)


def current_coord(xyz):
    theview = cmd.get_view()
    TMatrix = theview[:9]
    new = []
    for i in range(3):
        vector = [TMatrix[i], TMatrix[i+3], TMatrix[i+6]]
        tmp = sum([(xyz[ii] * v) for ii, v in enumerate(vector)])
        new.append(tmp)
    return new

def sele2CA_list(sele):
    my = cmd.get_model(sele)
    xyz = []
    for atom in my.atom:
        if atom.name == "CA":
            xyz.append(current_coord(atom.coord))
    return xyz

def rotateView(front, back):
    #Rotate 180 degree if necessary to put domain to the front
    atoms = sele2CA_list(front)
    num = len(atoms)
    if num == 0:
        print 'ERROR: No atoms found in domain'
        return False
    fz = sum([i[2] for i in atoms]) / num

    atoms = sele2CA_list(back)
    num = len(atoms)
    bz = sum([i[2] for i in atoms]) / num
    if fz < bz:
        cmd.turn('y', 180)
    return True


uid, pdb, r = args[:3]
ligand = None
asb_range = []
if len(args) > 3:
    options = args[3:]
    for opt in options:
        if not '-' in opt:
            #ligand annotation supposed to be single indexes.
            if ligand is not None:
                print 'WARNING: muliple ligand annotations given?  Overwriting...'
            ligand = opt
        else:
            asb_range.append(opt)


################################
#Load structures and preparation

cmd.set('ignore_case', '0')
cmd.set('retain_order', '1')
cmd.set('cartoon_smooth_loops', 0)
cmd.set('cartoon_side_chain_helper', 1)
cmd.set('sphere_scale', 0.5)
cmd.bg_color('white')

d = os.path.join(basepath, uid[2:-2], uid)
if not os.path.exists(d):
    print 'WARNING: domain dir does not exist: %s' % uid
    os.makedirs(d)
os.chdir(d)


pdbfile = pdbpath % (pdb[1:3], pdb)
if not os.path.exists(pdbfile):
    pdbfile = cifpath % (pdb[1:3], pdb)
    if not os.path.exists(pdbfile):
        print 'ERROR: pdb %s does not exists' % pdb
        sys.exit()


sele, all_chains = range2sele(r)
if pdbfile.endswith('.cif.gz'):
    import tempfile
    f = gzip.open(pdbfile, 'rb')
    tmpfile, tmpfilepath = tempfile.mkstemp(".cif", "tmpcif")
    f_out = os.fdopen(tmpfile, 'w')
    f_out.write(f.read())
    f_out.close()
    f.close()
    cmd.load(tmpfilepath, format='cif', object=pdb)
    os.remove(tmpfilepath)
else:
    #pdb.gz does not need decompress
    cmd.load(pdbfile, format='pdb', object=pdb)

cmd.color('grey80', 'all')
cmd.hide('everything', 'all')
cmd.select('dom', sele)
# Only multichain domain use custom rainbow function
if len(all_chains) > 1:
    rainbow(r)
else:
    cmd.spectrum('count', 'rainbow', 'dom')



#################################
#Show lignads in sticks or spheres
sphere_sele = []
stick_sele = []
if ligand is not None:
    lig_sele = ''
    pymol.stored.name = ''
    for lig in ligand.split(','):
        c, i = lig.split(':')
        if i.startswith('-'):
            # minus index need be escaped
            i = '\\' + i
        lig_sele = 'c. %s & i. %s' % (c, i)
        cmd.iterate(lig_sele, 'stored.name=resn')
        if pymol.stored.name in atomic_ligands:
            sphere_sele.append(lig_sele)
        else:
            stick_sele.append(lig_sele)
    if sphere_sele:
        sphere_sele = ' | '.join(sphere_sele)
        cmd.show('spheres', sphere_sele)
        #util.cbag(sphere_sele)
        cmd.color('atomic', sphere_sele)
        #show coordination side chains
        cmd.select('inter', 'byres (!(n. C+O+H+CA+N)&dom within 4 of (%s))' % sphere_sele)
        cmd.color('atomic','!e. C & inter & (!(n. C+O+H+CA|(n. N&!r. pro)))')
        cmd.color('grey70','e. C & inter & (!(n. C+O+H+CA|(n. N&!r. pro)))')
        cmd.show('sticks', 'inter & !(n;C,O,H|(n. N&!r. pro))')
    if stick_sele:
        stick_sele = ' | '.join(stick_sele)
        cmd.show('sticks', stick_sele)
        #color carbon with grey70
        cmd.color('atomic','!e. C & (%s)' % stick_sele)
        cmd.color('grey70','e. C & (%s)' % stick_sele)


#################################
# Image for domain in pdb context

cmd.select('ca', 'all and name CA')
ca_num = cmd.count_atoms('ca')
n_num = cmd.count_atoms('name N and bto. ca')
if n_num * 1.0 / ca_num < 0.33:
    cmd.set('ribbon_trace_atoms', 1)
    rep = 'ribbon'
else:
    rep = 'cartoon'
natoms = cmd.count_atoms('all')
if natoms > 99999:
    rep = 'ribbon'
if natoms <= 499999:
    cmd.orient('all')
    # zoom completely ensures all atoms in the view
    cmd.zoom('all', complete=1)
    cmd.show(rep, 'all')
    # show disulfide bonds, see pymol/menu.py
    cmd.show('sticks', '(byres ((dom & r. CYS+CYX & n. SG) &\
              bound_to (dom & r. CYS+CYX & n. SG))) & n.\
              CA+CB+SG')
    cmd.deselect()
    cmd.png('%s_pdb_thumb.png' % uid, 300, 300, ray=1)
    cmd.png('%s_pdb.png' % uid, 1024, 768, ray=1)
else:
    #special solution for large structures
    cmd.orient('dom')
    # zoom completely ensures all atoms in the view
    cmd.zoom('all', complete=1)
    cmd.color('red', 'dom')
    rotateView('dom', pdb)
    cmd.set('gaussian_resolution', 10)
    cmd.set('gaussian_b_floor', 50)
    cmd.map_new('map', 'gaussian', 10, pdb, 10)
    cmd.isosurface('surf', 'map', 1)
    cmd.ramp_new('ramp', pdb, [0,10,10], [-1,-1,0])
    cmd.color('ramp', 'surf')
    cmd.disable('ramp')
    cmd.set('transparency', 0.5)
    cmd.png('%s_pdb_thumb.png' % uid, 300, 300, ray=1)
    cmd.png('%s_pdb.png' % uid, 1024, 768, ray=1)
    cmd.delete('surf')
    if len(all_chains) > 1:
        rainbow(r)
    else:
        cmd.spectrum('count', 'rainbow', 'dom')
    cmd.set('transparency', 1)


##################################
# Image for domain in chain context

cmd.hide(rep, 'all')
cmd.set('ribbon_trace_atoms', 0)
chain_sele = ' | '.join(['c. %s' % i for i in all_chains])
cmd.select('ca', '(%s) and name CA' % chain_sele)
ca_num = cmd.count_atoms('ca')
if ca_num == 0:
    print 'WARNING: No CA atoms found!'
    rep = 'ribbon'
else:
    n_num = cmd.count_atoms('name N and bto. ca')
    if n_num * 1.0 / ca_num < 0.33:
        cmd.set('ribbon_trace_atoms', 1)
        rep = 'ribbon'
    else:
        rep = 'cartoon'
if cmd.count_atoms(chain_sele) > 99999:
    # although it should only exceed in whole PDB
    rep = 'ribbon'
cmd.show(rep, chain_sele)
cmd.orient(chain_sele)
comb_sele = chain_sele
if stick_sele:
    comb_sele += ' | %s' % stick_sele
if sphere_sele:
    comb_sele += ' | %s' % sphere_sele
cmd.zoom(comb_sele, complete=1)
cmd.deselect()
cmd.png('%s_chain_thumb.png' % uid, 300, 300, ray=1)
cmd.png('%s_chain.png' % uid, 1024, 768, ray=1)



######################
# Image for the domain

cmd.hide(rep, 'all')
cmd.set('ribbon_trace_atoms', 0)
if asb_range:
    # domain sele could change here, so domain image is the last
    sele += ''.join([' | %s' % range2sele(i)[0] for i in asb_range])
    cmd.save('%s_assembly.pdb' % uid, sele)
cmd.select('ca', '(%s) and name CA' % sele)
ca_num = cmd.count_atoms('ca')
n_num = cmd.count_atoms('name N and bto. ca')
if n_num * 1.0 / ca_num < 0.33:
    cmd.set('ribbon_trace_atoms', 1)
    rep = 'ribbon'
else:
    rep = 'cartoon'
if cmd.count_atoms(sele) > 99999:
    rep ='ribbon'
cmd.show(rep, sele)
cmd.orient(sele)
comb_sele = sele
if stick_sele:
    comb_sele += ' | %s' % stick_sele
if sphere_sele:
    comb_sele += ' | %s' % sphere_sele
cmd.zoom(comb_sele, complete=1)
cmd.deselect()
cmd.png('%s_thumb.png' % uid, 300, 300, ray=1)
cmd.png('%s.png' % uid, 1024, 768, ray=1)


# vim: ts=4 expandtab sw=4 sts=4 tw=78
