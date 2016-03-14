# This file is part of BratUtils.
#
# BratUtils is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BratUtils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BratUtils.  If not, see <http://www.gnu.org/licenses/>.
__author__ = 'Aleksandar Savkov'

"""This module groups data structure classes necessary for the calculation of
inter-annotator agreement (IAA) of two sets of parallel annotations. The main
statics counting class is MucTable, which the rest provide suitable brat enabled
data structures for it.
"""

import os
import glob
import ntpath


class Comparison:
    """A container for the parameters used in the comparison of two annotations.
    """

    def __init__(self):
        """All attributes are initialised as `False`.


        """
        self.borders = False
        self.tag = False
        self.partial = False

    @property
    def is_incorrect(self):
        """True if the `tag` attribute, or both `border` and `partial`
        attributes are set to `False`.


        :return: `result of comparison
        :rtype: bool
        """
        if (self.borders or self.partial) and self.tag:
            return False
        return True

    def __str__(self):
        return str({"borders": self.borders,
                    "tag": self.tag,
                    "partial": self.partial})


class MucTable:
    """Data structure object reflecting the parameters of the MUC-7 scoring
    scheme.
    """

    RELAXED_COMPARISON = 1
    STRICT_COMPARISON = 2
    BORDER_COMPARISON = 3

    CORRECT = 1
    INCORRECT = 2
    PARTIAL = 3
    MISSING = 4
    SPURIOUS = 5

    _tsvkeys = ["pos", "act", "cor", "par", "inc", "mis", "spu", "pre", "rec",
                "fsc", "und", "ovg", "sub", "bor", "ibo"]

    def __init__(self):
        self.comparison = None

        #MUC-7 scores

        self.pos = 0
        self.act = 0
        self.cor = 0
        self.par = 0
        self.inc = 0
        self.mis = 0
        self.spu = 0
        self.rec = 0.0
        self.pre = 0.0
        self.und = 0.0
        self.ovg = 0.0
        self.sub = 0.0
        self.fsc = 0.0

        self.bor = 0  # border correct
        self.ibo = 0  # incorrect border

    def update_table(self, comparison_type=None):
        """Updates the values of atributes of this object (e.g. `fsc`) that are
        calculated based on the values of the counted attributes, such as `cor`
        (correct).

        :param comparison_type: comparison type
        :type: one of STRICT_COMPARISON, RELAXED_COMPARISON or BORDER_COMPARISON
        """
        self.pos = self.cor + self.inc + self.mis
        self.act = self.cor + self.inc + self.spu
        self.bor = self.cor + self.inc
        self.ibo = self.par + self.spu
        if comparison_type is None:
            comparison_type = MucTable.STRICT_COMPARISON
        self.comparison = comparison_type
        if comparison_type == self.RELAXED_COMPARISON:
            self.pos = self.cor + self.par + self.inc + self.mis
            self.act = self.cor + self.par + self.inc + self.spu
            try:
                self.rec = float(self.cor + self.par) / self.pos
            except ZeroDivisionError:
                self.rec = 0.0
            try:
                self.pre = float(self.cor + self.par) / self.act
            except ZeroDivisionError:
                self.pre = 0.0
            try:
                self.und = float(self.mis) / self.pos
            except ZeroDivisionError:
                self.und = 0.0
            try:
                self.ovg = float(self.spu) / self.act
            except ZeroDivisionError:
                self.ovg = 0.0
            try:
                self.sub = float(self.inc) / (self.cor + self.par + self.inc)
            except ZeroDivisionError:
                self.sub = 0.0
        elif comparison_type == self.STRICT_COMPARISON:
            self.pos = self.cor + self.par + self.inc + self.mis
            self.act = self.cor + self.par + self.inc + self.spu
            try:
                self.rec = float(self.cor) / self.pos
            except ZeroDivisionError:
                self.rec = 0.0
            try:
                self.pre = float(self.cor) / self.act
            except ZeroDivisionError:
                self.pre = 0.0
            try:
                self.und = float(self.mis) / self.pos
            except ZeroDivisionError:
                self.und = 0.0
            try:
                self.ovg = float(self.spu) / self.act
            except ZeroDivisionError:
                self.ovg = 0.0
            try:
                self.sub = (float((self.inc + self.par)) /
                           (self.cor + self.par + self.inc))
            except ZeroDivisionError:
                self.sub = 0.0
        # TODO: add exception handling
        elif comparison_type == self.BORDER_COMPARISON:
            self.pos = self.bor + self.par + self.mis
            self.act = self.bor + self.par + self.spu
            self.rec = float(self.bor) / self.pos
            self.pre = float(self.bor) / self.act
            self.und = float(self.mis) / self.pos
            self.ovg = float(self.spu) / self.act
            self.sub = float(self.ibo) / (self.bor + self.ibo)
        else:
            print("Something's wrong!")

        try:
            self.fsc = 2.0 * (float(float(self.pre) * float(self.rec)) /
                              float(float(self.pre) + float(self.rec)))
        except ZeroDivisionError:
            self.fsc = 0.0

    def add_table(self, muc_table):
        """Accumulates the values of the counted attributes of another MucTable
        object into this one and updates the values of calculated attributes.

        :param muc_table: MucTable object
        :type: MucTable
        """
        self.pos += muc_table.pos
        self.act += muc_table.act
        self.cor += muc_table.cor
        self.par += muc_table.par
        self.inc += muc_table.inc
        self.mis += muc_table.mis
        self.spu += muc_table.spu
        self.bor += muc_table.bor
        self.ibo += muc_table.ibo
        self.update_table(None)

    @property
    def tsvheader(self):
        """The list of MUC-7 attributes separated with tabs. Should be used as
        header line in a TSV representation of this object.


        :return: tsv header
        :rtype: str
        """
        return "\t".join(self._tsvkeys)

    @property
    def tsvstring(self):
        """Tab separated value representation of the attributes of this object.


        :return: tsv string representation of this object
        :rtype: str
        """
        values = [str(self.__dict__.get(key)) for key in self._tsvkeys]
        return "\t".join(values)

    def __str__(self):
        keychains = [
            ["pos", "act", "cor", "par", "inc", "mis", "spu"],
            ["pre", "rec", "fsc"],
            ["und", "ovg", "sub"],
            ["bor", "ibo"]
        ]
        blocks = []
        line = "\n------------------------------------------------\n"
        title_line = "\n-------------------MUC-Table--------------------"
        blocks.append(title_line)
        for keychain in keychains:
            blocks.append(
                "\n".join(["%s:%s" % (str(key), str(self.__dict__.get(key)))
                           for key
                           in keychain]))
        blocks.append(line)
        return line.join(blocks)

    def __repr__(self):
        return str(self)


class Annotation:
    """Annotation data structure encoding the tag and position of an annotaiton,
    along with information about its comparison status.
    """

    def __init__(self, a):
        """Constructs an Annotation object from a brat annotation line string.

        :param a: annotation line string
        :type: str
        """
        self.text = None
        self.start_idx = None
        self.end_idx = None
        self.tag_name = None
        self.partial_match = None

        self.comp_status = None
        self.comp_match = None
        self.border_status = False
        self.border_match = None

        self.text, self.tag_name, self.start_idx, self.end_idx = \
            self._parse_annotation(a)

    @staticmethod
    def _parse_annotation(a):
        items = a.split("\t")
        text = items[2].strip("\n").strip(" ")
        subitems = items[1].split(" ")
        tag_name = subitems[0]
        start_idx = int(subitems[1])
        end_idx = int(subitems[2])
        return text, tag_name, start_idx, end_idx

    def reset_markers(self):
        """Resets the comparison marker attributes to default values. The
        default value for `comp_status` attribute is `SPURIOUS` in case it is
        not treated as gold standard.


        """
        self.partial_match = None
        self.comp_status = MucTable.SPURIOUS
        self.comp_match = None
        self.border_status = False
        self.border_match = None

    def make_gold(self):
        """Sets the `comp_status` attribute to MISSING, which is its default
        value when this object is considered gold standard.


        """
        self.comp_status = MucTable.MISSING

    def reverse_gold(self):
        """Reverses gold standard status and sets `comp_status` value again to
        `SPURIOUS`.


        """
        self.comp_status = MucTable.SPURIOUS

    def update_comp_status(self, comp):
        """Updates the comparison status attribute using a Comparison object.

        :param comp: comaprison
        :type comp: Comparison
        """
        if comp.borders:
            self.border_status = True
        else:
            self.border_status = False
        if comp.borders and comp.tag:
            self.comp_status = MucTable.CORRECT
        elif comp.partial:
            self.comp_status = MucTable.PARTIAL
        elif not comp.tag or not comp.borders:
            self.comp_status = MucTable.INCORRECT

    def compare_to(self, parallel_ann):
        """Compared this object to a parallel annotation.

        :param parallel_ann:
        :return: True if objects are the same
        :rtype: bool
        """
        comp = Comparison()
        if (self.start_idx == parallel_ann.start_idx and
                self.end_idx == parallel_ann.end_idx):
            comp.borders = True
        if self.tag_name == parallel_ann.tag_name:
            comp.tag = True
        if comp.borders and comp.tag:
            return comp
        tag_contained_in_gold = \
            (parallel_ann.start_idx <= self.start_idx < parallel_ann.end_idx or
             parallel_ann.start_idx < self.end_idx <= parallel_ann.end_idx)
        tag_span_over_gold = (self.start_idx < parallel_ann.start_idx and
                              self.end_idx > parallel_ann.end_idx)
        if ((tag_contained_in_gold or tag_span_over_gold) and
                self.tag_name == parallel_ann.tag_name and
                not parallel_ann.partiallyMatched):
            parallel_ann.partially_matched = True
            comp.partial = True
        return comp

    def coincides_with(self, parallel_ann):
        """Checks if this object's annotation coincides with the parallel
        annotation in another object.

        :param parallel_ann: parallel annotation object
        :return: True if objects coincide
        :rtype: bool
        """
        return (self.start_idx == parallel_ann.start_idx and
                self.end_idx == parallel_ann.end_idx)

    def contains_ann(self, other_ann):
        """Checks if this object's annotation contains another object's
        annotation.

        :param other_ann: annotation object
        :return: True if this annotaion contains the other annotation
        :rtype: bool
        """
        return (other_ann.start_idx >= self.start_idx and
                other_ann.end_idx <= self.end_idx)

    def is_contained_by(self, parallel_ann):
        """Checks if this annotation is contained by a parallel annotation.

        :param parallel_ann:
        :return: True if contained in `parallel_ann`
        :rtype: bool
        """
        return (parallel_ann.start_idx <= self.start_idx and
                parallel_ann.end_idx >= self.end_idx)

    def get_contained_anns(self, parallel_anns):
        """Returns a list of parallel annotations contained in this annotation.

        :param parallel_anns: list of parallel annotations
        :return: list of contained annotations
        :rtype: list
        """
        return [x for x in parallel_anns if self.contains_ann(x)]

    def get_containing_ann(self, parallel_anns):
        """Returns the parallel annotation that contains this annotation.

        :param parallel_anns: list of parallel annotations
        :return: containing parallel annotation
        :rtype: Annotation
        """
        it = iter(parallel_anns)
        try:
            t = it.next()
            while not self.is_contained_by(t):
                t = it.next()
        except StopIteration:
            return None
        return t

    def overlaps_with(self, parallel_ann):
        """Returns True if this annotation overlaps with the parallel annotaiton
        in `parallel_ann`.

        :param parallel_ann: parallel annotation
        :return: True if the annotations overlap
        :rtype: bool
        """
        if (parallel_ann.end_idx < self.start_idx or
                self.end_idx < parallel_ann.start_idx):
            return False
        elif self == parallel_ann:
            return True
        else:
            # TODO check indexes in cases like 'I'
            tag_contained_in_gold = (parallel_ann.start_idx <=
                                     self.start_idx <
                                     parallel_ann.end_idx or
                                     parallel_ann.start_idx <
                                     self.end_idx <=
                                     parallel_ann.end_idx)
            tag_span_over_gold = (self.start_idx <= parallel_ann.start_idx and
                                  self.end_idx >= parallel_ann.end_idx)
            return tag_contained_in_gold or tag_span_over_gold

    def get_overlapping_anns(self, parallel_anns):
        """Returns a list of parallel annotations overlapping with this one.

        :param parallel_anns: list of parallel annotations
        :return: overlapping annotations
        :rtype: list
        """
        overlapping_tags = []
        for gTag in parallel_anns:
            if self.overlaps_with(gTag):
                overlapping_tags.append(gTag)
        return overlapping_tags

    def has_partial_candidate(self, parallel_ann):
        """Checks the provided parallel annotation is a partial match candidate.

        :param parallel_ann: parallel annotation
        :return: True if parallel annotation is a partial match
        :rtype: bool
        """
        tag_contained_in_gold = (parallel_ann.start_idx <= self.start_idx and
                                 self.end_idx <= parallel_ann.end_idx)
        tag_span_over_gold = (self.start_idx <= parallel_ann.start_idx and
                              parallel_ann.end_idx <= self.end_idx)
        return tag_contained_in_gold or tag_span_over_gold

    def is_right_from(self, ann):
        """Checks if this annotation is right from another annotation `ann`.

        :param ann: another annotation
        :return:True if right from `ann`
        :rtype: bool
        """
        return ann.end_idx < self.start_idx

    def in_range(self, idx_range):
        """Checks if this annotation is in a specific index range `idx_range`.

        :param idx_range: index range
        :return: True if in the index range
        :rtype: bool
        """
        return (idx_range[0] <= self.start_idx <= idx_range[1] or
                idx_range[0] <= self.end_idx <= idx_range[1])

    def __eq__(self, ann):
        """Checks if this object is the same as `ann`. Objects are considered
        the same when they have the same text, start_idx, end_idx, and tag_name
        attributes.

        :param ann: another Annotation object
        :return: True if objects are the same
        :rtype: bool
        """
        return (self.text == ann.text and
                self.start_idx == ann.start_idx and
                self.end_idx == ann.end_idx and
                self.tag_name == ann.tag_name)

    def __str__(self):
        return " ".join([self.tag_name, self.text])

    def __repr__(self):
        return " ".join([self.tag_name, self.text])


class Filter:
    """Filters the entries in a Document object based on tag, borders or id.
    """

    TAG_FILTER = 'tag'
    BORDER_FILTER = 'border'
    ID_FILTER = 'id'

    def __init__(self, name, filter_type, scope, positive_polarity):

        self._filter_funcs = {
            self.TAG_FILTER: self._filter_tags,
            self.BORDER_FILTER: self._filter_borders,
            self.ID_FILTER: self._filter_ids
        }

        try:
            self.filter_func = self._filter_funcs[filter_type]
        except KeyError:
            print("Undefined filter type: %s" % filter_type)

        self.name = name
        self.conditions = scope
        self.positive_polarity = positive_polarity

    def apply_filter(self, document):
        """Applies this filter to `document`.

        :param document: Document object
        """
        self.filter_func(document)

    def _filter_tags(self, document):
        new_tags = []
        for tag in document.postag_list:
            for filter_tag in self.conditions:
                if ((self.positive_polarity and tag.tag_name == filter_tag) or
                        (not self.positive_polarity and
                            tag.tag_name is not filter_tag)):
                    new_tags.append(tag)
        document.postag_list = new_tags

    def _filter_borders(self, document):
        new_tags = []
        for tag in document.postag_list:
            for condition in self.conditions:
                in_range = tag.in_range(condition)
                if ((self.positive_polarity and not in_range) or
                        (not self.positive_polarity and in_range)):
                    new_tags.append(tag)
                else:
                    pass
        document.postag_list = new_tags

    def _filter_ids(self, document):
        pass


class Document:
    """A collection of Annotation objects usually originating from the same
    annotation document (.ann file).
    """

    def __init__(self, fp=None, ann_list=None):
        """Constructs a `Document` object from an annotation file using `fp` or
        using a collection of Annotation objects in `ann_list`.

        :param fp: annotation document file path
        :param ann_list: list of annotations
        """
        self.tags = []
        self.correct = []
        self.incorrect = []
        self.partial = []
        self.spurious = []
        self.missing = []
        if fp:
            with open(fp) as doc:
                for line in doc:
                    if not line.startswith("#"):
                        self.tags.append(Annotation(line))
        elif ann_list:
            for line in ann_list:
                if not line.startswith("#"):
                    self.tags.append(Annotation(line))
        else:
            self.tags = []
        self.sort()

    def sort(self):
        """Sort annotations in this document by their starting index.


        """
        self.tags.sort(key=lambda tag: tag.start_idx)

    def make_gold(self):
        """Set all annotations in this document to gold standard default values.
        Look at same method in Annotation.


        """
        for tag in self.tags:
            tag.make_gold()

    def reverse_gold(self):
        """Set all annotations in this document back to normal default values.
        Look at same method in Annotation.


        """
        for tag in self.tags:
            tag.reverse_gold()

    def compare_to_gold(self, parallel_doc):
        """Compares this annotation document to a parallel document.

        :param parallel_doc: parallel document
        :return: MucTable with degree of agreement
        :rtype: MucTable
        """
        muc = MucTable()
        self.sort()
        parallel_doc.sort()
        for tag in self.tags:
            contained_tags = tag.get_contained_anns(parallel_doc.tags)
            if not contained_tags:
                ctag = tag.get_containing_ann(parallel_doc.tags)
                if ctag and ctag.comp_status == MucTable.MISSING:
                    muc.par += 1
                    tag.comp_status = MucTable.PARTIAL
                    ctag.comp_status = MucTable.PARTIAL
                    continue
            for ctag in contained_tags:
                if tag == ctag:  # full match
                    tag.comp_status = MucTable.CORRECT
                    ctag.comp_status = MucTable.CORRECT
                    muc.cor += 1
                elif tag.coincides_with(ctag):  # incorrect tag
                    tag.comp_status = MucTable.INCORRECT
                    ctag.comp_status = MucTable.INCORRECT
                    muc.inc += 1
                elif (tag.comp_status != MucTable.PARTIAL and
                        tag.tag_name == ctag.tag_name):  # partial match
                    tag.comp_status = MucTable.PARTIAL
                    ctag.comp_status = MucTable.PARTIAL
                    muc.par += 1
                else:  # incorrect
                    if tag.comp_status != MucTable.PARTIAL:
                        tag.comp_status = MucTable.INCORRECT
                    ctag.comp_status = MucTable.INCORRECT
                    muc.inc += 1
        # Spurious tags are tags that are not fully contained in any gold tag
        # spu = n of guess tags - n of contained tags
        muc.spu = len(self.tags) - len([x for x in self.tags if x.comp_status])
        # Missing tags are tags that are not fully contained in any guess tag
        # mis = n of gold tags - n of contained tags
        muc.mis = len([x
                       for x
                       in parallel_doc.tags
                       if x.comp_status == MucTable.MISSING])
        muc.update_table()
        return muc

    def filter_document(self, doc_filters):
        """Applies the list of filter to this annotation document.

        :param doc_filters: list of filters
        :type doc_filters: list of Filter objects
        """
        for doc_filter in doc_filters:
            doc_filter.apply_filter(self)

    def reset_markers(self):
        """Reset markers of all annotations in this document.


        """
        for tag in self.tags:
            tag.reset_markers()

    def __str__(self):
        return "".join(self.tags)

    def __repr__(self):
        return "".join(self.tags)


class DocumentCollection:
    """A collection of annotation documents.
    """
    def __init__(self, document_dir_path, ext='ann'):
        """All annotation document files (default: .ann) located in a folder or
        one of its subfolders are included in this object.

        :param document_dir_path: directory with annotation document files
        :raise ValueError: Empty or non-existant folder
        """
        if not os.path.isdir(document_dir_path):
            raise IOError("The %s does not exist." % document_dir_path)
        if os.listdir(document_dir_path) is []:
            raise ValueError("Empty Document Collection directory: %s"
                             % document_dir_path)
        self.documents = {}
        files = glob.glob(os.path.join(document_dir_path, "*.%s" % ext))
        for fileName in files:
            self.documents[ntpath.basename(fileName)] = \
                Document(fp=fileName)

    @property
    def correct(self):
        """List of lists of correct annoations from this document collection.


        :return: correct annotations
        """
        correct = []
        for doc in self.documents.values():
            correct.extend(doc.correct)
        return correct

    @property
    def incorrect(self):
        """List of lists of incorrect annoations from this document collection.


        :return: incorrect annotations
        """
        incorrect = []
        for doc in self.documents.values():
            incorrect.extend(doc.incorrect)
        return incorrect

    @property
    def partial(self):
        """List of lists of partial annoations from this document collection.


        :return: partial annotations
        """
        partial = []
        for doc in self.documents.values():
            partial.extend(doc.partial)
        return partial

    @property
    def spurious(self):
        """List of lists of spurious annoations from this document collection.


        :return: spurious annotations
        """
        spurious = []
        for doc in self.documents.values():
            spurious.extend(doc.spurious)
        return spurious

    @property
    def missing(self):
        """List of lists of missing annoations from this document collection.


        :return: missing annotations
        """
        missing = []
        for doc in self.documents.values():
            missing.extend(doc.missing)
        return missing

    def make_gold(self):
        """Set all annoations in all document as gold standard.


        """
        for doc in self.documents.values():
            doc.make_gold()

    def reverse_gold(self):
        """Reverse gold standard of all annoation in all document.


        """
        for doc in self.documents.values():
            doc.reverse_gold()

    def reset_markers(self):
        """Reset all markers of all annoations in all documents.


        """
        for doc in self.documents.values():
            doc.reset_markers()

    def compare_to_gold(self, gold_collection):
        """Compare this document collection to a parallel gold standard version
        of the same documents.

        :param gold_collection: gold standard parallel document collection
        :return: accumulated average agreement results
        :rtype: MucTable
        """
        muc = MucTable()
        gold_collection.make_gold()
        for key in self.documents.keys():
            doc = self.documents.get(key)
            gold = gold_collection.documents.get(key)
            if doc is not None and gold is not None:
                muc.add_table(doc.compare_to_gold(gold))
            else:
                print("Warning: document naming mismatch. " + key)
        gold_collection.reverse_gold()
        return muc

    def filter_document_collection(self, filters):
        """Apply a filter list to all documents in this collection.

        :param filters: filter list
        """
        for document in self.documents.values():
            document.filter_document(filters)

    def __str__(self):
        return "".join(["%s\n\n%s\n" % (key, self.documents.get(key))
                        for key
                        in self.documents.keys()])

    def __repr__(self):
        return str(self)


