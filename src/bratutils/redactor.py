# This file is part of bratutils.
#
# bratutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# bratutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bratutils.  If not, see <http://www.gnu.org/licenses/>.


__author__ = 'Aleksandar Savkov'


def redact_line(l):
    its = l.split('\t')
    toks = its[-1].split(' ')
    redacted = []
    for t in toks:
        redacted.append('~' * len(t))
    its[-1] = ' '.join(redacted)
    return '\t'.join(its)


def redact_ann(fp, op):
    lines = []
    with open(fp, 'r') as fh:
        for l in fh:
            lines.append(redact_line(l))
    redacted = '\n'.join(lines)

    with open(op, 'w') as fh:
        fh.write(redacted)


def redact_text(fp, op):
    t = open(fp, 'r').read()
    lines = []
    for l in t.split('\n'):
        toks = []
        for tt in l.split(' '):
            toks.append('~' * len(tt))
        lines.append(' '.join(toks))
    redacted = '\n'.join(lines)

    with open(op, 'w') as fh:
        fh.write(redacted)
