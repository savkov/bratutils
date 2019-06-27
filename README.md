bratutils
=========
[![CircleCI](https://circleci.com/gh/savkov/bratutils.svg?style=svg&circle-token=9a7bdcb066c87c45017fe2214c71f2e2f9672c94)](https://circleci.com/gh/savkov/bratutils)
[![Maintainability](https://api.codeclimate.com/v1/badges/4c8fccbe0c29026c90bd/maintainability)](https://codeclimate.com/github/savkov/bratutils/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/4c8fccbe0c29026c90bd/test_coverage)](https://codeclimate.com/github/savkov/bratutils/test_coverage)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A collection of utilities for manipulating data and calculating inter-annotator 
agreement in brat annotation files.

### Installation

Install as a normal package from the source directory.

```bash
$ pip install bratutils
```


### Agreement Definition

Agreement in multi-token annotations is commonly evaluated using [f-score][fsc].
due to various problems with computing the traditional [Krippendorf's alpha][al] 
and [Cohen's kappa][ka]. [Hripcsak][hripcsak] prove the validity of the metric 
for very large populations, i.e. for unrestricted text annotations.

This library roughly follows the definitions of precision and recall calculation
from the [MUC-7 test scoring][muc]. The basic definitions along with some 
additional restrictions are laid out below:

* `CORRECT` - when annotation tags and indices match completely
* `INCORRECT` - when annotation tags do not match, but the indices coincide
* `PARTIAL` - when the annotation tags are the same but one of the annotations
has the same end index and a different start index
* `MISSING` - annotations exising only in the gold standard annotation set
* `SPURIOUS` - annotations existing only in the candidate annotation set

_Note_: the gold standard is considered the collections/document from which the 
 comparison is invoked, while the supplied parallel annotation is considered 
 the candidate set.
 
_*Disclaimer:*_ the current definition of the `PARTIAL` category accomodates 
working with syntactic chunks. A different arrangement (e.g. pick largest 
contained tag as partial match instead of rightmost) might be more suitable for 
other tasks, for example some types of semantic annotation.


### Examples

Simple example:

```python
from bratutils import agreement as a

doc = a.Document('res/samples/A/data-sample-1.ann')
doc2 = a.Document('res/samples/B/data-sample-1.ann')

doc.make_gold()
statistics = doc2.compare_to_gold(doc)

print(statistics)
```

Output:

```shell
-------------------MUC-Table--------------------
------------------------------------------------
pos:135
act:134
cor:115
par:5
inc:4
mis:11
spu:10
------------------------------------------------
pre:0.858208955224
rec:0.851851851852
fsc:0.855018587361
------------------------------------------------
und:0.0814814814815
ovg:0.0746268656716
sub:0.0725806451613
------------------------------------------------
bor:119
ibo:15
------------------------------------------------
------------------------------------------------
```


[fsc]: <https://en.wikipedia.org/wiki/F1_score>
[al]: <https://en.wikipedia.org/wiki/Krippendorff%27s_alpha>
[ka]: <https://en.wikipedia.org/wiki/Cohen%27s_kappa>
[hripcsak]: <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1090460/>
[muc]: <https://aclweb.org/anthology/M/M98/M98-1024.pdf>
