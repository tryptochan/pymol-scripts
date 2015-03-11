import re
import colorsys
from pymol import cmd


def rainbow(range_string):
    """
DESCRIPTION

    Colors rainbow spectrum for a selection given in range string.

    The difference between coloring in rainbow with built-in 'spectrum' is that
    this relies on the segment order in range string (not alphabetically
    sorted), so it can handle multiple chain domain as in insulin where usually
    chain B should be before chain A in many cases.

USAGE

    rainbow range_string

ARGUMENTS

    range_string = 'B:2-29,A:1-21'
    """
    seg_ptn = re.compile(r'([A-Za-z0-9]{1}):(-?[0-9]+[A-Z]?)-(-?[0-9]+[A-Z]?)')
    all_resi = []
    for seg in seg_ptn.finditer(range_string):
        chain = seg.group(1)
        local_space = {'resnums' : [], 'chain': chain}
        groups = list(seg.groups())
        for i in [1, 2]:
            # excape minus index
            if groups[i].startswith('-'):
                groups[i] = '\\' + groups[i]
        cmd.iterate('c. %s and i. %s-%s and n. CA' % seg.groups(),
                    'resnums.append(resi)', space=local_space)
        all_resi.append(local_space)

    total = reduce(lambda x, y: x + len(y['resnums']), all_resi, 0)

    cnt = 0
    for seg in all_resi:
        chain = seg['chain']
        for i in seg['resnums']:
            hue = colorsys.TWO_THIRD - colorsys.TWO_THIRD * cnt / (total - 1)
            red, green, blue = colorsys.hsv_to_rgb(hue, 1, 1)
            hexcolor = hex((int(red * 255) << 16) + (int(green * 255) << 8) +
                           int(blue * 255))
            cmd.color(hexcolor, 'c. %s and i. %s' % (chain, i))
            cnt += 1

if __name__ != "rainbow":
    cmd.extend('rainbow', rainbow)
