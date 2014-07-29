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
import os
from os.path import dirname
import shutil
from bratutils.data import BratDocument
from bratutils.agreement import MucTable

_penn_escape_dict = {'DOT': '.', 'OQ': '``', 'CQ': "''", 'CLM': ':', 'CM': ',', 'LPE': '(', 'RPE': ')'}


# Merges the annotation files across all subdirectories of the a and b, deploying the results in c
def merge_brat_dox_dir(dir_path_a, dir_path_b, dir_path_c):  # TODO: integrate with dir processor

    logging.info("Merging directories: %s & %s", dir_path_a, dir_path_b)
    logging.info("Target path: %s", dir_path_c)

    text = []
    annotation_files = []

    for root, dirs, files in os.walk(dir_path_a):

        for f in files:

            if f.endswith("txt"):
                text.append("{0}{1}{2}".format(root, os.path.sep, f))

            elif f.endswith("ann"):
                annotation_files.append("{0}{1}{2}".format(root, os.path.sep, f))

    for ann_file_path in set(annotation_files):

        ann_dir = dir_path_a if dir_path_a in ann_file_path else dir_path_b
        mirrored_ann_file_path = ann_file_path.replace(dir_path_a, dir_path_b) if dir_path_a in ann_file_path else ann_file_path.replace(dir_path_b, dir_path_a)
        merged_ann_file_path = ann_file_path.replace(ann_dir, dir_path_c)
        text_file = "%stxt" % ann_file_path[0:-3]
        text_file_copy = "%stxt" % merged_ann_file_path[0:-3]

        if not os.path.exists(dirname(merged_ann_file_path)):
            os.makedirs(dirname(merged_ann_file_path))

        logging.debug("Text file: %s Text file copy: %s", text_file, text_file_copy)
        logging.debug("Annotation1: %s\nAnnotation2: %s\nMerged: %s\n", ann_file_path, mirrored_ann_file_path, merged_ann_file_path)

        merge_brat_documents(ann_file_path, mirrored_ann_file_path, merged_ann_file_path)
        shutil.copyfile(text_file, text_file_copy)

    logging.info("Merging done.")


# Merges two annotation of the same document (a and b) into their union without duplicates stored into document c.
def merge_brat_documents(doc_path_a, doc_path_b, doc_path_c):
    doc1 = BratDocument(doc_path=doc_path_a)
    doc2 = BratDocument(doc_path=doc_path_b)
    merged_doc = BratDocument()
    merged_doc.update(doc1)
    doc1size = len(doc1)
    doc1comments = doc1.get_num_comments()
    for record in doc2.values():
        record_id = int(record.id)
        record.update_id(record_id + doc1size)
        if record.comment:
            record.comment.id = "{0}{1}".format(record.comment.id, doc1comments)
        merged_doc[record.id] = record
    merged_doc.remove_duplicates()
    merged_doc.enumerate_comments()
    merged_doc.export_to_file(doc_path_c)


def unescape_brat_tags(ann_file_path, escape_dict=_penn_escape_dict):
    ann = BratDocument(doc_path=ann_file_path)
    ann.unescape_tags(escape_dict)
    ann.export_to_file(ann_file_path)


class MissingAnnotationFileError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def get_fscore_vector(muc_tables, param="fsc", comparison=MucTable.HARD_COMPARISON):
    vector = []
    for muc in muc_tables:
        muc.update_table(comparison)
        vector.append(muc.__dict__[param])
    return vector