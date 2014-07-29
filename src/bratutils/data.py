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

import logging


class BratComment(object):

    def __init__(self, commentstr):
        items = commentstr.split("\t")
        self.id = int(items[0][1:])
        self.note = items[2]
        self.recordref = int(items[1][16:])

    def __str__(self):
        return "#{0}\tAnnotatorNotes T{1}\t{2}\n".format(str(self.id),
                                                         str(self.recordref),
                                                         self.note)

    def __eq__(self, other):
        if not isinstance(other, BratComment):
            return False
        if self.id != other.id:
            return False
        if self.note != other.note:
            return False
        if self.recordref != other.recordref:
            return False
        return True

    def __hash__(self):
        return hash((self.id, self.recordref, self.note))


class BratAnnotation(object):
    def __init__(self, recordstr):
        items = recordstr.split("\t")
        self.id = int(items[0][1:])
        subitems = items[1].split(" ", 1)
        self.tag = subitems[0]
        self.boundaries = []
        for boundarystr in subitems[1].split(";"):
            self.boundaries.append(boundarystr.split(" "))
        self.content = items[2]
        self.comment = None

    def update_id(self, new_id):
        self.id = new_id
        if self.comment:
            self.comment.recordref = new_id

    def set_comment(self, comment):
        self.comment = comment

    def add_comment(self, comment):
        if self.comment and comment:
            self.comment.note = "Comment 1: {0} Comment 2: {1}".format(
                self.comment.note.replace("\n", ""), comment.note)
        elif comment:
            self.comment = comment
            self.comment.recordref = self.id

    def get_boundary_text(self):
        return ";".join([" ".join(x) for x in self.boundaries])

    def get_left_border(self):
        return int(self.boundaries[0][0])

    def get_right_border(self, inclusive=False):
        return (int(self.boundaries[-1][1])-1
                if inclusive
                else int(self.boundaries[-1][1]))

    def __str__(self):
        comment_str = str(self.comment) if self.comment else ""
        return "T{0}\t{1} {2}\t{3}\n{4}".format(str(self.id), self.tag,
                                                ";".join([" ".join(b)
                                                          for b
                                                          in self.boundaries]),
                                                self.content, comment_str)

    def __repr__(self):
        comment_str = "; {0}".format(str(self.comment)) if self.comment else ""
        return "T{0}\t{1} {2}\t{3}{4}".format(str(self.id),
                                              self.tag,
                                              ";".join([" ".join(b)
                                                        for b
                                                        in self.boundaries]),
                                              self.content, comment_str)

    def __eq__(self, other):
        if not isinstance(other, BratAnnotation):
            return False
        if self.id != other.id:
            return False
        if self.tag != other.tag:
            return False
        if self.content != other.content:
            return False
        if not self.comment == other.comment:
            return False
        if self.boundaries != other.boundaries:
            return False
        return True

    def __hash__(self):
        return hash((self.id, self.boundaries, self.content, self.comment))


class BratDocument(dict):

    def __init__(self, doc_path=None, no_newline=True):
        super(BratDocument, self).__init__()
        if doc_path:
            with open(doc_path) as doc:
                for line in doc:
                    if no_newline:
                        line = line[0:-1]
                    if line.startswith("#"):
                        comment = BratComment(line)
                        self[comment.recordref].set_comment(comment)
                    else:
                        ann = BratAnnotation(line)
                        self[ann.id] = ann

    def filter_tags(self, filters, positive_polarity=True):
        for key in self.keys():
            if (self[key].tag in filters) != positive_polarity:
                del self[key]

    def unescape_tags(self, escaped_tags):
        for ann in self.values():
            if ann.tag in escaped_tags.keys():
                ann.tag = escaped_tags[ann.tag]

    def remove_duplicates(self):
        duplicates = []
        for key in self.keys():
            if key in duplicates:
                continue
            new_duplicates = self.get_duplicates(key)
            for duplicate in new_duplicates:
                self[key].add_comment(self[duplicate].comment)
                del self[duplicate]
            duplicates.extend(new_duplicates)

    def get_duplicates(self, target_key):
        duplicates = []
        for key in self.keys():
            if key == target_key:
                continue
            if self[key].boundaries == (self[target_key].boundaries and
                                        self[key].tag == self[target_key].tag):
                duplicates.append(key)
        return duplicates

    def enumerate_comments(self):
        comment_id = 1
        for record in self.values():
            if record.comment:
                record.comment.id = comment_id
                comment_id += 1

    def get_num_comments(self):
        count = 0
        for record in self.values():
            if record.comment:
                count += 1
        return count

    def get_sorted_values(self):
        return sorted(self.values(), cmp=_comp_brat_record)

    def validate_indices(self):
        triples = []
        for record in self.values():
            triples.append((record.get_left_border(), True, record))
            triples.append((record.get_right_border(), False, record))
        triples.sort()
        for i, tr in enumerate(triples):
            look_ahead = 1
            balance = 0
            stack = []
            while tr[1] and tr[2] != triples[i + look_ahead][2]:
                if triples[i + look_ahead][1]:
                    balance += 1
                else:
                    balance -= 1
                stack.append(triples[i + look_ahead])
                look_ahead += 1
            if balance != 0:
                print "Stack %s\n%s" % (tr, stack)
                logging.info("Stack %s\n%s" % (tr, stack))
                return False
        return True

    def export_to_file(self, path):
        with open(path, "w") as doc:
            for record in self.values():
                doc.write(str(record))

    def __eq__(self, other):
        if not isinstance(other, BratDocument):
            return False
        if self.keys() != other.keys():
            return False
        for key in self.keys():
            if not self[key] == other[key]:
                return False
        return True


def _comp_brat_record(record_a, record_b):

    hash_a = (record_a.get_left_border(), record_a.get_right_border())
    hash_b = (record_b.get_left_border(), record_b.get_right_border())

    if hash_a < hash_b:
        return 1
    elif hash_a == hash_b:
        return 0
    else:
        return -1

