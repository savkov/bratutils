import os
import shutil
import logging

from os.path import dirname
from bratutils.bratdata import BratAnnotation
from bratutils.agreement import MucTable


__author__ = 'Aleksandar Savkov'

"""Utility package for the bratutils package.
"""


_penn_escape_dict = {'DOT': '.', 'OQ': '``', 'CQ': "''", 'CLM': ':', 'CM': ',',
                     'LPE': '(', 'RPE': ')'}


class MissingAnnotationFileError(Exception):
    """Thrown when one of the brat annotation file pair is missing. Each brat
     annotation consists of a text `txt` file and an annotation `ann` file of
     the same name (before the file extension suffix).
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def merge_brat_dox_dir(dp_a, dp_b, dp_res):
    """Merges the annotation files across all subdirectories of `dp_a` and
    `dp_b`, deploying the results in `dp_res`.

    :param dp_a: directory path A
    :param dp_b: directory path B
    :param dp_res: results directory path
    """
    logging.info("Merging directories: %s & %s", dp_a, dp_b)
    logging.info("Target path: %s", dp_res)

    text = []
    annotation_files = []

    for root, dirs, files in os.walk(dp_a):

        for f in files:

            if f.endswith("txt"):
                text.append("{0}{1}{2}".format(root, os.path.sep, f))

            elif f.endswith("ann"):
                annotation_files.append(
                    "{0}{1}{2}".format(root, os.path.sep, f))

    for ann_file_path in set(annotation_files):

        ann_dir = dp_a if dp_a in ann_file_path else dp_b
        mirrored_ann_file_path = (ann_file_path.replace(dp_a,
                                                        dp_b)
                                  if dp_a in ann_file_path
                                  else ann_file_path.replace(dp_b,
                                                             dp_a))
        merged_ann_file_path = ann_file_path.replace(ann_dir, dp_res)
        text_file = "%stxt" % ann_file_path[0:-3]
        text_file_copy = "%stxt" % merged_ann_file_path[0:-3]

        if not os.path.exists(dirname(merged_ann_file_path)):
            os.makedirs(dirname(merged_ann_file_path))

        logging.debug("Text file: %s Text file copy: %s",
                      text_file,
                      text_file_copy)
        logging.debug("Annotation1: %s\nAnnotation2: %s\nMerged: %s\n",
                      ann_file_path,
                      mirrored_ann_file_path,
                      merged_ann_file_path)

        merge_brat_documents(ann_file_path,
                             mirrored_ann_file_path,
                             merged_ann_file_path)
        shutil.copyfile(text_file, text_file_copy)

    logging.info("Merging done.")


def merge_brat_documents(fp_a, fp_b, fp_res):
    """Merges two annotation of the same document (`fp_a` amd `fp_b`) into
    their union without duplicates stored into document in `fp_res`.

    :param fp_a: file path A
    :param fp_b: file path B
    :param fp_res: results file path
    """
    # TODO: this looks wrong; investigate
    doc1 = BratAnnotation(doc_path=fp_a)
    doc2 = BratAnnotation(doc_path=fp_b)
    merged_doc = BratAnnotation()
    merged_doc.update(doc1)
    doc1_size = len(doc1)
    doc1_comments = doc1.comments_count
    for record in doc2.values():
        record_id = int(record.id)
        record.update_id(record_id + doc1_size)
        if record.comment:
            record.comment.id = "{0}{1}".format(record.comment.id,
                                                doc1_comments)
        merged_doc[record.id] = record
    merged_doc.remove_duplicates()
    merged_doc.enumerate_comments()
    merged_doc.export_to_file(fp_res)


def get_stat_vector(muc_tables,
                    param="fsc",
                    comparison=MucTable.STRICT_COMPARISON):
    """Returns a statistic vector extracted from a list of MucTable objects.

    Note: this method re-estimates the statistics in the `MucTable` object.

    :param muc_tables: list of `MucTable` objects
    :param param: statistic name
    :param comparison: comparison type
    :return: a statistic vector
    :rtype: list
    """
    vector = []
    for muc in muc_tables:
        muc.update_table(comparison)
        vector.append(muc.__dict__[param])
    return vector
