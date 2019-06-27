__author__ = 'Aleksandar Savkov'

"""This module holds data structure classes and relevant methods for parsing
and manipulating brat annotation files.0
"""

import logging


class BratComment(object):
    """Brat comment container
    """

    def __init__(self, commentstr):
        """Creates a new object from a comment line from a brat annotation file.

        :param commentstr: comment line string
        """

        self.id, self.note, self.recordref = self._parse_comment(commentstr)

    @staticmethod
    def _parse_comment(c):
        """Parses comment line string into id, note, and reference

        :param c: comment line string
        :return: id, note, and reference
        :rtype: str, str, str
        """
        items = c.split("\t")
        rid = items[0]
        note = items[2]
        ref = items[1].rsplit(' ', 1)[-1]
        return rid, note, ref

    def __str__(self):
        return "{0}\tAnnotatorNotes {1}\t{2}\n".format(str(self.id),
                                                       str(self.recordref),
                                                       self.note)

    def __eq__(self, other):
        """Returns `True` if the other object is of type `BratComment` and the
        values of the three attributes `id`, `recordref`, and `note` are
        identical.

        :param other:
        :return: True if the objects are the same
        :rtype: bool
        """
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
        """Returns hash based on the `id`, `recordref`, and `note` attributes.


        :return: hash
        """
        return hash((self.id, self.recordref, self.note))


class BratAnnotation(object):
    """Brat annotation usually recorded on a single line in a .ann file.
    """

    def __init__(self, s):
        self.comment = None
        self.id, self.tag, self.boundaries, self.content = \
            self._parse_annotation(s)

    @staticmethod
    def _parse_annotation(a):
        items = a.split("\t")
        rid = items[0]
        subitems = items[1].split(" ", 1)
        tag = subitems[0]
        boundaries = []
        for boundarystr in subitems[1].split(";"):
            boundaries.append(boundarystr.split(" "))
        content = items[2]
        return rid, tag, boundaries, content

    def update_id(self, new_id):
        """Update a record ID and any attached to it comments

        :param new_id: new id
        """
        self.id = new_id
        if self.comment:
            self.comment.recordref = new_id

    def set_comment(self, comment):
        """Set comment value

        :param comment: comment
        """
        self.comment = comment

    def add_comment(self, comment):
        """Add new comment to this annotation

        :param comment:
        """
        if self.comment and comment:
            self.comment.note = "Comment 1: {0} Comment 2: {1}".format(
                self.comment.note.replace("\n", ""), comment.note)
        elif comment:
            self.comment = comment
            self.comment.recordref = self.id

    @property
    def boundary_str(self):
        """Returns string representation of the boundaries of this annoation.


        :return:
        """
        return ";".join([" ".join(x) for x in self.boundaries])

    def get_left_border(self):
        """Returns left border index of this annotation.


        :return: left border index
        :rtype: int
        """
        return int(self.boundaries[0][0])

    def get_right_border(self, inclusive=False):
        """Returns right border index of this annotation.

        :param inclusive: inclusive border position
        :return: right border index
        :rtype: int
        """
        return (int(self.boundaries[-1][1])-1
                if inclusive
                else int(self.boundaries[-1][1]))

    def __str__(self):
        comment_str = str(self.comment) if self.comment else ""
        return "{0}\t{1} {2}\t{3}\n{4}".format(str(self.id), self.tag,
                                               ";".join([" ".join(b)
                                                         for b
                                                         in self.boundaries]),
                                               self.content, comment_str)

    def __repr__(self):
        comment_str = "; {0}".format(str(self.comment)) if self.comment else ""
        return "{0}\t{1} {2}\t{3}{4}".format(str(self.id),
                                             self.tag,
                                             ";".join([" ".join(b)
                                                       for b
                                                       in self.boundaries]),
                                             self.content, comment_str)

    def __eq__(self, other):
        """Returns True if the following attributes of both objects are the
        same: `id`, `tag`, `content`, `comment`, `boundaries`

        :param other: other annotation object
        :return: True if objects are the same
        :rtype: bool
        """
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
        """Returns a hash based on the following attributes: `id`, `tag`,
        `content`, `comment`, `boundaries`


        :return: hash code
        """
        return hash((self.id, self.boundaries, self.content, self.comment))


class BratDocument(dict):
    """Annotation file document made up of annotations represented with
    Annotation objects. This object extends dict. Annotations are stored using
    their ID as a key value.
    """

    def __init__(self, fp=None, no_newline=True, string=None):
        """Construct a new object and optionally add content to it from an
        annotation document file (.ann).

        :param fp: annotation document file path
        :param no_newline: trim newline char
        """
        super(BratDocument, self).__init__()
        if fp:
            with open(fp) as doc:
                self._parse_document(doc, no_newline=no_newline)
        elif string:
            self._parse_document(string.split('\n'), no_newline=no_newline)

    def _parse_document(self, doc, no_newline):
        """
        Parse an annotation document by iterating over its lines
        """
        for line in doc:
            if no_newline:
                line = line.strip()
            if line.startswith("#"):
                comment = BratComment(line)
                self[comment.recordref].set_comment(comment)
            elif line[0] in 'RAE':
                # annotations that are not handled ATM
                continue
            elif line.startswith('T'):
                ann = BratAnnotation(line)
                self[ann.id] = ann
            else:
                raise ValueError('Unknown beginning of line')

    def filter_tags(self, filters, positive_polarity=True):
        """Filter annotations in this document

        :param filters: filter list
        :param positive_polarity: keep filtered values
        """
        for key in self.keys():
            if (self[key].tag in filters) != positive_polarity:
                del self[key]

    def unescape_tags(self, escaped_tags_dict):
        """Replace escaped annotation tags with original values. Brat can't
        handle certain characters, such as dollar sign in the tag annotation.

        Dictionary example:
        {'NN_DOLLAR': 'NN$'}

        :param escaped_tags_dict: dictionary of escaped tags
        """
        for ann in self.values():
            if ann.tag in escaped_tags_dict.keys():
                ann.tag = escaped_tags_dict[ann.tag]

    def remove_duplicates(self):
        """Removes duplicate annotations in this document.


        """
        duplicates = []
        for key in self.keys():
            if key in duplicates:
                continue
            new_duplicates = self.get_duplicates(key)
            for duplicate in new_duplicates:
                self[key].add_comment(self[duplicate].comment)
            duplicates.extend(new_duplicates)
        for duplicate in duplicates:
            del self[duplicate]

    def get_duplicates(self, target_key):
        """Returns a list of keys of duplicate annotations of a target
        annotation.

        :param target_key: target key
        :return: keys of duplicate annotations
        :rtype: list
        """
        duplicates = []
        for key in self.keys():
            if key == target_key:
                continue
            if (self[key].boundaries == self[target_key].boundaries and
                    self[key].tag == self[target_key].tag):
                duplicates.append(key)
        return duplicates

    def enumerate_comments(self):
        """Enumerate all comments attached to annotations in this document.


        """
        comment_id = 1
        for record in self.values():
            if record.comment:
                record.comment.id = comment_id
                comment_id += 1

    @property
    def comments_count(self):
        """Number of comments in this document.


        :return: number of comments
        :rtype: int
        """
        count = 0
        for record in self.values():
            if record.comment:
                count += 1
        return count

    @property
    def sorted_values(self):
        """Returns a list of this documents annotations sorted by their border
        indices.


        :return: list of annotations
        :rtype: list
        """
        return sorted(self.values(), cmp=_comp_brat_ann)

    def validate_indices(self):
        """Validates the annotation indices in this document.


        :return: True if indices are correct
        :rtype: bool
        """
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
                print("Stack %s\n%s" % (tr, stack))
                logging.info("Stack %s\n%s" % (tr, stack))
                return False
        return True

    def export_to_file(self, fp):
        """Exports this document to a file.

        :param fp: file path
        """
        with open(fp, "w") as doc:
            for record in self.values():
                doc.write(str(record))

    def __eq__(self, other):
        """Returns True if the key sets and their respective values are the
        same.

        :param other:
        :return:
        """
        if not isinstance(other, BratDocument):
            return False
        if self.keys() != other.keys():
            return False
        for key in self.keys():
            if not self[key] == other[key]:
                return False
        return True


def _comp_brat_ann(ann_a, ann_b):
    """Compare brat annotation records by their border indices.

    :param ann_a: annotation A
    :param ann_b: annotation B
    :return: 1 a < b; 0 a == b; -1 a > b
    :rtype: int
    """
    hash_a = (ann_a.get_left_border(), ann_a.get_right_border())
    hash_b = (ann_b.get_left_border(), ann_b.get_right_border())

    if hash_a < hash_b:
        return 1
    elif hash_a == hash_b:
        return 0
    else:
        return -1
