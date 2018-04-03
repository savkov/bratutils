# This file is part of bratutils.
#
# bratutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# bratutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bratutils.  If not, see <http://www.gnu.org/licenses/>.
from bratutils import agreement as a


__author__ = 'Aleksandar Savkov'

# doc = a.Document('../res/samples/other/final_embedding_a.ann')
# doc2 = a.Document('../res/samples/other/final_embedding_b.ann')
#
# doc.make_gold()
# statistics = doc2.compare_to_gold(doc)
#
# print statistics

doc = a.Document('../res/samples/A/data-sample-1.ann')
doc2 = a.Document('../res/samples/B/data-sample-1.ann')

# doc.make_gold()
statistics = doc2.compare_to_gold(doc)

print statistics