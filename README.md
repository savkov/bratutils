BratUtils
=========

A collection of utilities for manipulating data and calculating inter-annotator 
agreement in brat annotation files.

### Installation

Install as a normal package from the source directory.

```python setup install```


### Examples

Simple example:

```python
from bratutils import agreement as a

doc = a.Document('res/samples/A/data-sample-1.ann')
doc2 = a.Document('res/samples/B/data-sample-1.ann')

doc.make_gold()
statistics = doc2.compare_to_gold(doc)

print statistics

-------------------MUC-Table--------------------

\------------------------------------------------

pos:158

act:158

cor:130

par:0

inc:28

mis:0

spu:0

\------------------------------------------------

pre:0.822784810127

rec:0.822784810127

fsc:0.822784810127

\------------------------------------------------

und:0.0

ovg:0.0

sub:0.177215189873

\------------------------------------------------

bor:158

ibo:0

\------------------------------------------------

\------------------------------------------------
```

### TODO

* migrate to Python 3
* add examples for `bratdata`
* document `utils`