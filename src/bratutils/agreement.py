import os
import glob
import ntpath
import logging


__author__ = 'Aleksandar Savkov'

"""This module groups data structure classes necessary for the calculation of
inter-annotator agreement (IAA) of two sets of parallel annotations. The main
statics counting class is MucTable, which the rest provide suitable brat
enabled data structures for it.
"""


def safe_division(a, b):
    """

    :param a: quantity A
    :param b: quantity B
    :return: the quatient
    :rtype: float
    """

    try:
        return a / b
    except ZeroDivisionError:
        return 0.0


def standard_logger(name='bratutils', log_path=None, log_level=logging.INFO):

    # create logger
    logger = logging.getLogger(name)

    # set level
    logger.setLevel(log_level)

    # create handlers
    h = logging.StreamHandler()
    h.setLevel(log_level)

    # set the format
    formatter = logging.Formatter('%(asctime)s :: %(message)s')
    h.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(h)

    # repeat for a file log
    if log_path:
        try:
            os.makedirs(log_path)
        except OSError:
            pass
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


logger = standard_logger(name='agreement', log_level=logging.DEBUG)


class Comparison:
    """A container for the parameters used in the comparison of two
    annotations.
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

        # MUC-7 scores

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
        :type: one of STRICT_COMPARISON, RELAXED_COMPARISON or
        BORDER_COMPARISON
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
            self.rec = safe_division(float(self.cor + self.par), self.pos)
            self.pre = safe_division(float(self.cor + self.par), self.act)
            self.und = safe_division(float(self.mis), self.pos)
            self.ovg = safe_division(float(self.spu), self.act)
            self.sub = safe_division(float(self.inc),
                                     (self.cor + self.par + self.inc))
        elif comparison_type == self.STRICT_COMPARISON:
            self.pos = self.cor + self.par + self.inc + self.mis
            self.act = self.cor + self.par + self.inc + self.spu
            self.rec = safe_division(float(self.cor), self.pos)
            self.pre = safe_division(float(self.cor), self.act)
            self.und = safe_division(float(self.mis), self.pos)
            self.ovg = safe_division(float(self.spu), self.act)
            try:
                self.sub = (float((self.inc + self.par)) /
                            (self.cor + self.par + self.inc))
            except ZeroDivisionError:
                self.sub = 0.0
        elif comparison_type == self.BORDER_COMPARISON:
            self.pos = self.bor + self.par + self.mis
            self.act = self.bor + self.par + self.spu
            self.rec = safe_division(float(self.bor), self.pos)
            self.pre = safe_division(float(self.bor), self.act)
            self.und = safe_division(float(self.mis), self.pos)
            self.ovg = safe_division(float(self.spu), self.act)
            self.sub = safe_division(float(self.ibo), (self.bor + self.ibo))
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
        self.update_table()

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
    """Annotation data structure encoding the tag and position of an
    annotation, along with information about its comparison status.
    """

    def __init__(self, a):
        """Constructs an Annotation object from a brat annotation line string.

        :param a: annotation line string
        :type: str
        """
        self.text = None
        self.frag = None            # fragment, usually one pair of indexes, but annotation can be
                                    # composed of several fragments (aka discontinuated annotations)
        self.start_idx = None       # start index of first fragment
        self.end_idx = None         # end index of last fragment
        self.tag_name = None
        self.partial_match = None

        self.comp_status = None
        self.comp_match = None
        self.border_status = False
        self.border_match = None

        self.text, self.tag_name, self.start_idx, self.end_idx, self.frag = \
            self._parse_annotation(a)

    @staticmethod
    def _parse_annotation(a):
        items = a.split("\t")
        text = items[2].strip("\n").strip(" ")
        subitems = items[1].split(" ")
        tag_name = subitems.pop(0)
        subitems = " ".join(subitems).split(";")
        frag = []
        for idx in subitems:
            start_idx, end_idx = idx.split(" ")
            frag.append((int(start_idx), int(end_idx)))
        start_idx = frag[0][0]
        end_idx = frag[len(frag)-1][1]
        return text, tag_name, start_idx, end_idx, frag

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

        return self.frag == parallel_ann.frag


    def contains_ann(self, other_ann):
        """Checks if this object's annotation contains another object's
        annotation.

        :param other_ann: annotation object
        :return: True if this annotaion contains the other annotation
        :rtype: bool
        """
        contained_fragments = [False] * len(other_ann.frag)
        for i in range(len(other_ann.frag)):
            for j in range(len(self.frag)):
                if other_ann.frag[i][0] >= self.frag[j][0] and \
                        other_ann.frag[i][1] <= self.frag[j][1]:
                    contained_fragments[i] = True
        # return True if all fragments of other_ann are contained in self
        return contained_fragments == [True] * len(other_ann.frag)


    def is_contained_by(self, parallel_ann):
        """Checks if this annotation is contained by a parallel annotation.

        :param parallel_ann:
        :return: True if contained in `parallel_ann`
        :rtype: bool
        """
        return self.contains_ann(parallel_ann)


    def is_partial_to(self, parallel_ann):
        """Returns `True` if the annotation is a partial match to the parallel
        annotation. To be considered a partial match the annotations have to
        have the same end index and tag name. The relation is not reflexive,
        therefore only the smaller annotation (smaller span) is considered to
        fulfill it.

        :param parallel_ann:
        :return:
        """
        # TODO really dive into frag (for now, we check start of first fragment and end of last fragment)
        return (self.start_idx > parallel_ann.end_idx and
                self.start_idx == parallel_ann.end_idx and
                self.tag_name == parallel_ann.tag_name)

    def get_same_anns(self, parallel_anns):
        """Returns a list of parallel annotations with that match this
        annotation.

        :param parallel_anns: list of parallel annotations
        :return: list of contained annotations
        :rtype: list
        """
        same = []
        for ann in parallel_anns:
            if self == ann:
                same.append(ann)
                logger.debug('Matching annotations: {} : {}'.format(self, ann))
        return same

    def get_coinciding_anns(self, parallel_anns):
        """Returns a list of annotations from a parallel annotation that
        conincide with the current annotation.

        :param parallel_anns: parallel annotation
        :return: list of coinciding annotations
        """
        coinciding = []
        for ann in parallel_anns:
            if self.coincides_with(ann):
                coinciding.append(ann)
                logger.debug('Coinciding annotations: {} : {}'
                             .format(self, ann))
        return coinciding

    def get_contained_anns(self, parallel_anns):
        """Returns a list of parallel annotations contained in this annotation.

        :param parallel_anns: list of parallel annotations
        :return: list of contained annotations
        :rtype: list
        """
        contained = []
        for ann in parallel_anns:
            if self.contains_ann(ann):
                contained.append(ann)
                logger.debug('`{}` contained in `{}`'.format(ann, self))
        return contained

    def get_containing_ann(self, parallel_anns):
        """Returns the parallel annotation that contains this annotation.

        :param parallel_anns: list of parallel annotations
        :return: containing parallel annotation
        :rtype: Annotation
        """

        containing_anns = [a for a in parallel_anns if self.is_contained_by(a)]
        containing_anns.sort(key=(lambda x: x.end_idx - x.start_idx))

        return containing_anns

    def overlaps_with(self, parallel_ann):
        """Returns True if this annotation overlaps with the parallel
        annotation in `parallel_ann`.

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
        """Checks the provided parallel annotation is a partial match
        candidate.

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
                self.frag == ann.frag and
                self.tag_name == ann.tag_name)

    def __str__(self):
        atts = [self.tag_name, str(self.start_idx), str(self.end_idx), str(self.frag),
                self.text]
        return " ".join(atts)

    def __repr__(self):
        atts = [self.tag_name, str(self.start_idx), str(self.end_idx), str(self.frag),
                self.text]
        return " ".join(atts)


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
        self.basename = ""
        if fp:
            self.basename = os.path.basename(fp)
            with open(fp, encoding='utf-8') as doc:
                for line in doc:
                    if not line.startswith("#") and not line.startswith("A"):  # ignoring Attributes
                        self.tags.append(Annotation(line))
        elif ann_list:
            for line in ann_list:
                if not line.startswith("#") and not line.startswith("A"):
                    self.tags.append(Annotation(line))
        else:
            self.tags = []
        self.sort()

    def sort(self):
        """Sort annotations in this document by their starting index.
        """

        self.tags.sort(key=lambda tag: tag.start_idx)

    def make_gold(self):
        """Set all annotations in this document to gold standard default
        values. Look at same method in Annotation.


        """
        for tag in self.tags:
            tag.make_gold()

    def reverse_gold(self):
        """Set all annotations in this document back to normal default values.
        Look at same method in Annotation.


        """
        for tag in self.tags:
            tag.reverse_gold()

    @staticmethod
    def handle_coinciding_tags(tag, ctags, muc):
        for ctag in ctags:
            if tag == ctag:
                tag.comp_status = MucTable.CORRECT
                ctag.comp_status = MucTable.CORRECT
                muc.cor += 1
                logger.debug('Correct match: {} : {}'.format(tag, ctag))
            else:
                ctag.comp_status = MucTable.INCORRECT
                muc.inc += 1
                logger.debug('Incorrect match: {} : {}'.format(tag, ctag))
        if tag.comp_status != MucTable.CORRECT:
            tag.comp_status = MucTable.INCORRECT

    @staticmethod
    def handle_contained_tags(tag, ctags, muc):
        for ctag in ctags:
            if ctag.is_partial_to(tag):
                par = 0
                if tag.comp_status != MucTable.CORRECT:
                    tag.comp_status = MucTable.PARTIAL
                    par = 1
                if ctag.comp_status != MucTable.CORRECT:
                    ctag.comp_status = MucTable.PARTIAL
                    par = 1
                muc.par += par
                logger.debug('Patrtial match: {} : {}'.format(tag, ctag))

    @staticmethod
    def handle_containing_tags(tag, ctags, muc):
        for ctag in ctags:
            logger.debug('`{}` :contains: `{}`'.format(ctag, tag))
            if tag.is_partial_to(ctag):
                par = 0
                if tag.comp_status != MucTable.CORRECT:
                    tag.comp_status = MucTable.PARTIAL
                    par = 1
                if ctag.comp_status != MucTable.CORRECT:
                    ctag.comp_status = MucTable.PARTIAL
                    par = 1
                muc.par += par
                logger.debug('Partial match: `{}` in `{}`'.format(tag, ctag))

    @staticmethod
    def count_remaining(parallel_annotations):
        remaining = []
        for ann in parallel_annotations:
            if ann.comp_status is None:
                remaining.append(ann)
        return len(remaining)

    def compare_to_gold(self, parallel_doc):
        """Compares this annotation document to a parallel document.

        :param parallel_doc: parallel document
        :return: MucTable with degree of agreement
        :rtype: MucTable
        """
        muc = MucTable()

        logger.debug('File: {}'.format(self.basename))
        logger.debug('Annotations A: {}\tAnnotations B: {}'
                     .format(len(self.tags), len(parallel_doc.tags)))

        self.sort()
        parallel_doc.sort()
        self.remove_duplicates()
        parallel_doc.remove_duplicates()
        for tag in self.tags:

            # tags with coinciding indices
            coinciding_tags = tag.get_coinciding_anns(parallel_doc.tags)
            if coinciding_tags:
                self.handle_coinciding_tags(tag, coinciding_tags, muc)
                continue

            # tags contained in the current annotation
            contained_tags = tag.get_contained_anns(parallel_doc.tags)
            if contained_tags:
                self.handle_contained_tags(tag, contained_tags, muc)
                continue

            # tags containing the current tag
            containing_tags = tag.get_containing_ann(parallel_doc.tags)
            if containing_tags:
                self.handle_containing_tags(tag, containing_tags, muc)
        # Spurious tags are tags that are not fully contained in any gold tag
        muc.spu = self.count_remaining(parallel_doc.tags)

        # Missing tags are tags that are not fully contained in any guess tag
        muc.mis = self.count_remaining(self.tags)

        muc.update_table()
        return muc

    def remove_duplicates(self):
        """Removes duplicate annotations in this document.
        """

        for tag in self.tags:
            equal_tags = tag.get_same_anns(self.tags)

            if len(equal_tags) > 1:
                for equal_tag in equal_tags[1:]:
                    self.tags.remove(equal_tag)

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
