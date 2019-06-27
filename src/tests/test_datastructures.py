from unittest import TestCase
from bratutils.bratdata import *

__author__ = 'Aleksandar Savkov'


class TestDataStructures(TestCase):

    def setUp(self):
        self.commentStr = "#1	AnnotatorNotes T36	Abbreviation for the " \
                          "word increase"
        self.comment2Str = "#2	AnnotatorNotes T45	Another abbreviation for" \
                           " the word increase"
        self.recordStr = "T36	MV 681 685	incr"
        self.record2Str = "T39	NP 798 805;810 819	bla bla"
        self.record3Str = "T45	NP 681 685	incr"

    def test_comment(self):
        c = BratComment(self.commentStr)
        self.assertEqual(c.id, '#1')
        self.assertEqual(c.recordref, 'T36')
        self.assertEqual(c.note, "Abbreviation for the word increase")

    def test_record(self):
        r = BratAnnotation(self.recordStr)
        self.assertEqual(r.id, 'T36')
        self.assertEqual(r.tag, "MV")
        self.assertEqual(r.boundaries, [["681", "685"]])
        self.assertEqual(r.content, "incr")

    def test_skip_relations(self):
        relation_ann = "E3	Positive_regulation:T4	Theme:E4"
        ann = BratDocument(string=relation_ann)
        self.assertEqual(len(ann), 0)

    def test_record_addcomment(self):
        r = BratAnnotation(self.recordStr)
        c = BratComment(self.commentStr)
        c2 = BratComment(self.comment2Str)
        r.set_comment(c)
        r.add_comment(c2)
        self.assertEqual(str(r.comment),
                          "#1\tAnnotatorNotes T36\tComment "
                          "1: Abbreviation for the word increase Comment "
                          "2: Another abbreviation for the word increase\n")

    def test_record_getBoundaryText(self):
        r = BratAnnotation(self.recordStr)
        r2 = BratAnnotation(self.record2Str)
        self.assertEqual(r.boundary_str, "681 685")
        self.assertEqual(r2.boundary_str, "798 805;810 819")

    def test_record_eq(self):
        r1 = BratAnnotation(self.recordStr)
        r2 = BratAnnotation(self.recordStr)
        self.assertEqual(r1, r2)
        self.assertListEqual([r1], [r2])

    # def test_tokenvector_getTokenIdx(self):
    #     tv = BratTokenVector(doc_path="./res/sampleTokenVector.tok")
    #     token = BratToken("symptomAs\t21\t29")
    #     idx = tv.get_token_idx(token)
    #     self.assertEqual(idx, 2)
    #
    # def test_tokenvector_insertListAt(self):
    #     tv = BratTokenVector(doc_path="./res/sampleTokenVector.tok")
    #     tokenList = [BratToken("foo	10	12"), BratToken("bar	21	23")]
    #     tok_before = BratToken("distension	10	19")
    #     tok_after = BratToken("symptomAs	21	29")
    #     tv.insert_list_at(2, tokenList)
    #     self.assertEqual(tv[1], tok_before, "Preceding token misatch.")
    #     self.assertEqual(tv[4], tok_after, "Following token misatch.")
    #     self.assertEqual(tv[2:4], tokenList, "List tokens mismatch.")

    def test_document_eq(self):
        doc = BratDocument()
        doc2 = BratDocument()
        r1 = BratAnnotation("T1	NP 0 28	Abdominal distension symptom\n")
        r11 = BratAnnotation("T1	NP 0 28	Abdominal distension symptom\n")
        c1 = BratComment(
            "#1	AnnotatorNotes T1	Abbreviation for the word increase #1\n")
        c11 = BratComment(
            "#1	AnnotatorNotes T1	Abbreviation for the word increase #1\n")
        r1.set_comment(c1)
        r11.set_comment(c11)
        doc[r1.id] = r1
        doc2[r11.id] = r11
        r2 = BratAnnotation("T2	TE 32 41	As before\n")
        r22 = BratAnnotation("T2	TE 32 41	As before\n")
        c2 = BratComment("#2	AnnotatorNotes T2	Another abbreviation for "
                         "the word increase #2\n")
        c22 = BratComment("#2	AnnotatorNotes T2	Another abbreviation for "
                          "the word increase #2\n")
        r2.set_comment(c2)
        r22.set_comment(c22)
        doc[r2.id] = r2
        doc2[r22.id] = r22
        self.assertEqual(r1, r11)
        self.assertEqual(r2, r22)
        self.assertEqual(doc, doc2)

    def test_document(self):
        doc_from_file = BratDocument("res/tests/sampledoc.ann")
        doc = BratDocument()
        r1 = BratAnnotation("T1	NP 0 28	Abdominal distension symptom")
        c1 = BratComment("#1	AnnotatorNotes T1	Abbreviation for the "
                         "word increase #1")
        r1.set_comment(c1)
        doc[r1.id] = r1
        r2 = BratAnnotation("T2	TE 32 41	As before")
        c2 = BratComment("#2	AnnotatorNotes T2	Another abbreviation "
                         "for the word increase #2")
        r2.set_comment(c2)
        doc[r2.id] = r2
        r3 = BratAnnotation("T3	NP 43 60	fybogel unhelpful")
        doc[r3.id] = r3
        r4 = BratAnnotation("T4	NP 0 28	Abdominal distension symptom")
        c3 = BratComment("#3	AnnotatorNotes T4	Another abbreviation for "
                         "the word increase #3")
        r4.set_comment(c3)
        doc[r4.id] = r4
        r5 = BratAnnotation("T5	TE 67 77	for a week")
        doc[r5.id] = r5
        r6 = BratAnnotation("T6	OE 79 86	On exam")
        doc[r6.id] = r6
        r7 = BratAnnotation("T7	NP 43 60	fybogel unhelpful")
        c4 = BratComment("#3	AnnotatorNotes T7	Another abbreviation for "
                         "the word increase #4")
        r7.set_comment(c4)
        doc[r7.id] = r7
        self.assertDictEqual(doc, doc_from_file)

    def test_document_get_duplicates(self):
        doc_from_file = BratDocument("res/tests/sampledoc.ann")
        duplicates = doc_from_file.get_duplicates('T1')
        self.assertListEqual(duplicates, ['T4'])
        duplicates = doc_from_file.get_duplicates('T2')
        self.assertListEqual(duplicates, [])
        duplicates = doc_from_file.get_duplicates('T3')
        self.assertListEqual(duplicates, ['T7'])
        duplicates = doc_from_file.get_duplicates('T4')
        self.assertListEqual(duplicates, ['T1'])

    def test_document_removeduplicates(self):
        doc_from_file = BratDocument("res/tests/sampledoc.ann")
        doc = BratDocument()
        r1 = BratAnnotation("T1	NP 0 28	Abdominal distension symptom")
        c1 = BratComment(
            "#1	AnnotatorNotes T1	"
            "Comment 1: Abbreviation for the word increase #1 "
            "Comment 2: Another abbreviation for the word increase #3")
        r1.set_comment(c1)
        doc[r1.id] = r1
        r2 = BratAnnotation("T2	TE 32 41	As before")
        c2 = BratComment("#2	AnnotatorNotes T2	Another abbreviation for "
                         "the word increase #2")
        r2.set_comment(c2)
        doc[r2.id] = r2
        r3 = BratAnnotation("T3	NP 43 60	fybogel unhelpful")
        c4 = BratComment("#3	AnnotatorNotes T3	Another abbreviation for "
                         "the word increase #4")
        r3.set_comment(c4)
        doc[r3.id] = r3
        r5 = BratAnnotation("T5	TE 67 77	for a week")
        doc[r5.id] = r5
        r6 = BratAnnotation("T6	OE 79 86	On exam")
        doc[r6.id] = r6
        doc_from_file.remove_duplicates()
        self.assertDictEqual(doc, doc_from_file)

    def test_unsupported_annotations(self):
        doc = BratDocument("res/tests/unsupporteddoc.ann")
        self.assertEqual(9, len(doc.values()))
