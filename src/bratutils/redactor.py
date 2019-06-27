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
