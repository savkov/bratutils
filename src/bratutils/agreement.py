# This file is part of DrTokenizer.
#
# DrTokenizer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DrTokenizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DrTokenizer.  If not, see <http://www.gnu.org/licenses/>.
__author__ = 'Aleksandar Savkov'

import glob
import os
import ntpath
import warnings


class Comparison:
    def __init__(self):
        self.borders = False
        self.tag = False
        self.partial = False

    def is_incorrect(self):
        if (self.borders or self.partial) and self.tag:
            return False
        return True

    def __str__(self):
        return str({"borders": self.borders,
                    "tag": self.tag,
                    "partial": self.partial})


class MucTable:
    SOFT_COMPARISON = 1
    HARD_COMPARISON = 2
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
        self.pos = self.cor + self.inc + self.mis
        self.act = self.cor + self.inc + self.spu
        self.bor = self.cor + self.inc
        self.ibo = self.par + self.spu
        if comparison_type is None:
            comparison_type = MucTable.HARD_COMPARISON
        self.comparison = comparison_type
        if comparison_type == self.SOFT_COMPARISON:
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
        elif comparison_type == self.HARD_COMPARISON:
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

    def print_out(self):
        warnings.warn("Method is deprecated. Use str() instead.",
                      DeprecationWarning)
        keys = ["pos", "act", "cor", "par", "inc", "mis", "spu"]
        keys2 = ["pre", "rec", "fsc"]
        keys3 = ["und", "ovg", "sub"]
        keys4 = ["bor", "ibo"]
        print("------------------------------------------------")
        for key in keys:
            print(str(key) + ":" + str(self.__dict__.get(key)))
        print("------------------------------------------------")
        for key in keys2:
            print(str(key) + ":" + str(self.__dict__.get(key)))
        print("------------------------------------------------")
        for key in keys3:
            print(str(key) + ":" + str(self.__dict__.get(key)))
        print("------------------------------------------------")
        for key in keys4:
            print(str(key) + ":" + str(self.__dict__.get(key)))
        print("------------------------------------------------\n")

    def get_tsvheader(self):
        return "\t".join(self._tsvkeys)

    def get_tsvstring(self):
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
        line = "------------------------------------------------"
        title_line = "-------------------MUC-Table--------------------"
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


class Tag:
    def __init__(self, tag_record):
        self.text = None
        self.start_idx = None
        self.end_idx = None
        self.tag_name = None
        self.partial_match = None

        self.comp_status = None
        self.comp_match = None
        self.border_status = False
        self.border_match = None

        tag_record_items = tag_record.split("\t")
        self.text = tag_record_items[2].strip("\n").strip(" ")
        tag_record_sub_items = tag_record_items[1].split(" ")
        self.tag_name = tag_record_sub_items[0]
        self.start_idx = int(tag_record_sub_items[1])
        self.end_idx = int(tag_record_sub_items[2])

    def reset_markers(self):
        self.partial_match = None
        self.comp_status = MucTable.SPURIOUS
        self.comp_match = None
        self.border_status = False
        self.border_match = None

    def make_gold(self):
        self.comp_status = MucTable.MISSING

    def reverse_gold(self):
        self.comp_status = MucTable.SPURIOUS

    def update_comp_status(self, comp):
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

    def compare_to(self, gold_tag):
        comp = Comparison()
        if (self.start_idx == gold_tag.start_idx and
                self.end_idx == gold_tag.end_idx):
            comp.borders = True
        if self.tag_name == gold_tag.tag_name:
            comp.tag = True
        if comp.borders and comp.tag:
            return comp
        tag_contained_in_gold = \
            (gold_tag.start_idx <= self.start_idx < gold_tag.end_idx or
             gold_tag.start_idx < self.end_idx <= gold_tag.end_idx)
        tag_span_over_gold = (self.start_idx < gold_tag.start_idx and
                              self.end_idx > gold_tag.end_idx)
        if ((tag_contained_in_gold or tag_span_over_gold) and
                self.tag_name == gold_tag.tag_name and
                not gold_tag.partiallyMatched):
            gold_tag.partially_matched = True
            comp.partial = True
        return comp

    def coincide_with(self, gold_tag):
        return (self.start_idx == gold_tag.start_idx and
                self.end_idx == gold_tag.end_idx)

    def contains_tag(self, gold_tag):
        return (gold_tag.start_idx >= self.start_idx and
                gold_tag.end_idx <= self.end_idx)

    def is_contained_by(self, gold_tag):
        return (gold_tag.start_idx <= self.start_idx and
                gold_tag.end_idx >= self.end_idx)

    def get_contained_tags(self, gold_tags):
        return [x for x in gold_tags if self.contains_tag(x)]

    def get_containing_tag(self, gold_tags):
        it = iter(gold_tags)
        try:
            t = it.next()
            while not self.is_contained_by(t):
                t = it.next()
        except StopIteration:
            return None
        return t

    def overlaps_with(self, gold_tag):
        if (gold_tag.end_idx < self.start_idx or
                self.end_idx < gold_tag.start_idx):
            return False
        elif self == gold_tag:
            return True
        else:
            # TODO check indexes in cases like 'I'
            tag_contained_in_gold = gold_tag.start_idx <= self.start_idx < gold_tag.end_idx or gold_tag.start_idx < self.end_idx <= gold_tag.end_idx
            tag_span_over_gold = self.start_idx <= gold_tag.start_idx and self.end_idx >= gold_tag.end_idx
            return tag_contained_in_gold or tag_span_over_gold

    def get_overlapping_tags(self, gold_tags):
        overlapping_tags = []
        for gTag in gold_tags:
            if self.overlaps_with(gTag):
                overlapping_tags.append(gTag)
        return overlapping_tags

    def has_partial_candidate(self, gold_tag):
        tag_contained_in_gold = gold_tag.start_idx <= self.start_idx and self.end_idx <= gold_tag.end_idx  # or goldTag.start_idx < self.end_idx <= goldTag.end_idx
        tag_span_over_gold = self.start_idx <= gold_tag.start_idx and gold_tag.end_idx <= self.end_idx
        return tag_contained_in_gold or tag_span_over_gold

    def is_right_from(self, tag):
        return tag.end_idx < self.start_idx

    def contained_in(self, tags):
        return self in tags

    def deep_contained_in(self, tags):
        for tag in tags:
            if self == tag:
                return True
        return False

    def in_range(self, idx_range):
        return idx_range[0] <= self.start_idx <= idx_range[1] or idx_range[0] <= self.end_idx <= idx_range[1]

    def __eq__(self, tag):
        return self.text == tag.text and self.start_idx == tag.start_idx and self.end_idx == tag.end_idx and self.tag_name == tag.tag_name

    def __str__(self):
        return " ".join([self.tag_name, self.text])

    def __repr__(self):
        return " ".join([self.tag_name, self.text])


class Filter:

    TAG_FILTER = 1
    BORDER_FILTER = 2
    ID_FILTER = 3

    def __init__(self, name, filter_type, scope, positive_polarity):
        self.name = name
        self.type = filter_type
        self.conditions = scope
        self.positive_polarity = positive_polarity

    def apply_filter(self, document):
        if self.type == self.TAG_FILTER or self.type == "tag":
            self.filter_tags(document)
        elif self.type == self.BORDER_FILTER or self.type == "border":
            self.filter_borders(document)
        elif self.type == self.ID_FILTER or self.type == "id":
            self.filter_ids(document)
        else:
            print("Undefined filter type!")

    def filter_tags(self, document):
        new_tags = []
        for tag in document.postag_list:
            for filter_tag in self.conditions:
                if (self.positive_polarity and tag.tag_name == filter_tag) or (not self.positive_polarity and tag.tag_name is not filter_tag):
                    new_tags.append(tag)
        document.postag_list = new_tags

    def filter_borders(self, document):
        new_tags = []
        for tag in document.postag_list:
            for condition in self.conditions:
                in_range = tag.in_range(condition)
                if (self.positive_polarity and not in_range) or (not self.positive_polarity and in_range):
                    new_tags.append(tag)
                else:
                    pass
        document.postag_list = new_tags

    def filter_ids(self, document):
        pass


class Document:  # the sum of all entries from an '.ann' file
    def __init__(self, document_path=None, document_list=None):
        self.tags = []
        self.correct = []
        self.incorrect = []
        self.partial = []
        self.spurious = []
        self.missing = []
        if document_path:
            with open(document_path) as doc:
                for line in doc:
                    if not line.startswith("#"):
                        self.tags.append(Tag(line))
        elif document_list:
            for line in document_list:
                if not line.startswith("#"):
                    self.tags.append(Tag(line))
        else:
            self.tags = []
        self.sort()

    def sort(self):
        self.tags.sort(key=lambda tag: tag.start_idx)

    def make_gold(self):
        for tag in self.tags:
            tag.make_gold()

    def reverse_gold(self):
        for tag in self.tags:
            tag.reverse_gold()

    def check_for_spurious(self, gold):
        mismatches = []
        for tag in self.tags:
            mismatch = False
            if tag.text == "able":
                pass
            for gTag in gold.postag_list:
                if gTag.start_idx <= tag.start_idx <= gTag.end_idx or gTag.start_idx <= tag.end_idx <= gTag.end_idx:
                    mismatch = False
                    break
                else:
                    mismatch = True
            if mismatch:
                mismatches.append(tag)
        return mismatches

    def check_for_missing(self, gold):
        return gold.check_for_spurious(self)

    def compare_to_gold(self, gold):
        '''
        Tags in this object are assumed to be guess tags. Gold standard tags are passed as parameter gold.
        '''
        muc = MucTable()
        self.sort()
        gold.sort()
        for tag in self.tags:
            contained_tags = tag.get_contained_tags(gold.postag_list)
            if not contained_tags:
                ctag = tag.get_containing_tag(gold.postag_list)
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
                elif tag.coincide_with(ctag):  # incorrect tag
                    tag.comp_status = MucTable.INCORRECT
                    ctag.comp_status = MucTable.INCORRECT
                    muc.inc += 1
                elif tag.comp_status != MucTable.PARTIAL and tag.tag_name == ctag.tag_name:  # partial match
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
        muc.mis = len([x for x in gold.postag_list if x.comp_status == MucTable.MISSING])
        muc.update_table()
        return muc

    def compare_to_gold_old(self, gold):
        muc = MucTable()
        for tag in self.tags:
            overlapping_tags = tag.get_overlapping_tags(gold.postag_list)
            partial_tag = None
            for overlapping_tag in overlapping_tags:
                if tag == overlapping_tag:  # the tag is a full match: cor
                    tag.comp_status = MucTable.CORRECT
                    overlapping_tag.comp_status = MucTable.CORRECT
                    self.correct.append(tag)
                    muc.cor += 1
                    break
                elif tag.coincide_with(overlapping_tag):  # the tag has matching borders, but wrong tag: inc
                    tag.comp_status = MucTable.INCORRECT
                    overlapping_tag.comp_status = MucTable.INCORRECT
                    self.incorrect.append(tag)
                    muc.inc += 1
                    break
                elif tag.tag_name == overlapping_tag.tag_name and tag.has_partial_candidate(overlapping_tag) and (partial_tag is None or overlapping_tag.is_right_from(partial_tag)):  # potential partial match: par
                    old_partial_match = overlapping_tag.partial_match
                    if old_partial_match is not None:
                        if tag.is_right_from(old_partial_match):
                            old_partial_match.comp_status = MucTable.SPURIOUS
                            if old_partial_match in self.partial:
                                self.partial.remove(old_partial_match)
                                muc.par -= 1
                                muc.spu += 1
                        else:
                            continue
                    partial_tag = overlapping_tag
                else:
                    # TODO: confirm if this should be muc.inc too
                    muc.inc += 1
            if partial_tag:
                muc.par += 1
                self.partial.append(tag)
                tag.comp_status = MucTable.PARTIAL
                partial_tag.comp_status = MucTable.PARTIAL
                partial_tag.partial_match = tag
                continue
            if tag.comp_status == MucTable.SPURIOUS:
                muc.spu += 1
                self.spurious.append(tag)
        for gold_tag in gold.postag_list:
            if gold_tag.comp_status == MucTable.MISSING:
                muc.mis += 1
        muc.update_table()
        return muc

    def filter_document(self, doc_filters):
        for doc_filter in doc_filters:
            doc_filter.apply_filter(self)

    def reset_markers(self):
        for tag in self.tags:
            tag.reset_markers()

    def __str__(self):
        return "".join(self.tags)

    def __repr__(self):
        return "".join(self.tags)


class DocumentCollection:  # collection of the entries from '.ann' files
    def __init__(self, document_dir_path):
        if not os.path.isdir(document_dir_path):
            raise IOError("The %s does not exist." % document_dir_path)
        if os.listdir(document_dir_path) is []:
            raise ValueError("Empty Document Collection directory: %s" % document_dir_path)
        self.documents = {}
        files = glob.glob(os.path.join(document_dir_path, "*.ann"))
        for fileName in files:
            self.documents[ntpath.basename(fileName)] = Document(document_path=fileName)

    def get_correct(self):
        correct = []
        for doc in self.documents.values():
            correct.extend(doc.correct)
        return correct

    def get_incorrect(self):
        incorrect = []
        for doc in self.documents.values():
            incorrect.extend(doc.incorrect)
        return incorrect

    def get_partial(self):
        partial = []
        for doc in self.documents.values():
            partial.extend(doc.partial)
        return partial

    def get_spurious(self):
        spurious = []
        for doc in self.documents.values():
            spurious.extend(doc.spurious)
        return spurious

    def get_missing(self):
        missing = []
        for doc in self.documents.values():
            missing.extend(doc.missing)
        return missing

    def make_gold(self):
        for doc in self.documents.values():
            doc.make_gold()

    def reverse_gold(self):
        for doc in self.documents.values():
            doc.reverse_gold()

    def reset_markers(self):
        for doc in self.documents.values():
            doc.reset_markers()

    def compare_to_gold(self, gold_collection):
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
        for document in self.documents.values():
            document.filter_document(filters)

    def to_token_vector_list(self):
        tv_list = []
        for document in self.documents.values():
            tv_list.append(document.toTokenVector())
        return tv_list

    def __str__(self):
        return "".join(["%s\n\n%s\n" % (key, self.documents.get(key)) for key in self.documents.keys()])

    def __repr__(self):
        return str(self)


