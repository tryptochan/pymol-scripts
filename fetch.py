from pymol import cmd
from pymol import preset
import re
import urllib2


def fetchd(code, name='', type='pdb', quiet=1, async=-1, **kwargs):
    """
DESCRIPTION

    Wraps and extends the built-in fetch function.

    It can fetch SCOP domains or PDBs and render the domain or chain of
    interest to publication reprentation.

USAGE

    fetch 1ast
    fetch 1astA
    fetch d1asta_
    """
    scop_url = 'http://scop.berkeley.edu/astral/pdbstyle/?id=%s&output=html'
    ptn_chain = re.compile(r"[0-9]{1}[A-Za-z0-9]{3}[A-Za-z0-9]{1}")
    ptn_scop = re.compile(r"d[a-z0-9]{4}[A-Za-z0-9.]{1}[0-9_]{1,2}")
    from pymol.importing import fetch as pymol_fetch

    quiet, async = int(quiet), int(async)
    if type != 'pdb':
        return pymol_fetch(code, name, type=type, quiet=quiet,
                           async=async, **kwargs)


    code_list = code.split()
    for code in list(code_list):
        if len(code) == 5 and ptn_chain.match(code):
            print 'Fetching PDB chain...'
            pdb_id = code[:4]
            chain_id = code[4]
            pymol_fetch(pdb_id, name, type=type, quiet=quiet, async=0, **kwargs)
            cmd.hide('everything')
            cmd.select('chain%s' % chain_id, '%s and chain %s' % (pdb_id, chain_id))
            cmd.deselect()
            preset.publication('chain%s' % chain_id)
        elif len(code) in (7, 8) and ptn_scop.match(code):
            print 'Fetching SCOP domain...'
            if 'proxy' in kwargs:
                proxy_handler = urllib2.ProxyHandler({'http': kwargs['proxy']})
                opener = urllib2.build_opener(proxy_handler)
            else:
                opener = urllib2.build_opener()
            content = opener.open(scop_url % code).read()
            cmd.read_pdbstr(content.replace('\n', '\\\n'), code)
            preset.publication("%s" % code)
        else:
            print "Normal pymol fetch..."
            return pymol_fetch(' '.join(code_list), name, quiet=quiet, async=async, **kwargs)

fetchd.__doc__ += cmd.fetch.__doc__
cmd.extend('fetch', fetchd)
