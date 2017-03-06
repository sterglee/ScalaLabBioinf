/*
 * Bioinformatics Algorithms in Scala
 * Copyright (C) 2016  Jason Mar
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bioinf.molecularevolution

import bioinf.UnitSpec
import MolecularEvolution._
import bioinf.mutations.Mutations.AdjacencyList

class MolecularEvolutionSpec extends UnitSpec {
  "MolecularEvolution" should "find distance between leaves" in {
    val v = IndexedSeq[Int](0, 1,2,3, 4,4,4,5,5,5)
    val w = IndexedSeq[Int](4, 4,5,5, 0,1,5,4,3,2)
    val l = IndexedSeq[Int](11,2,6,7,11,2,4,4,7,6)
    val x = AdjacencyList(v,w,l)
    val fw = floydWarshallAlgorithm(x)

    getDistanceBetweenLeaves(0,0,x,fw) should be (0)
    getDistanceBetweenLeaves(0,1,x,fw) should be (13)
    getDistanceBetweenLeaves(0,2,x,fw) should be (21)
    getDistanceBetweenLeaves(0,3,x,fw) should be (22)
    getDistanceBetweenLeaves(1,2,x,fw) should be (12)
    getDistanceBetweenLeaves(1,3,x,fw) should be (13)
    getDistanceBetweenLeaves(2,3,x,fw) should be (13)
  }

  it should "create additive phylogeny" in {
    val n = 4
    val D = DistanceMatrix(IndexedSeq[IndexedSeq[Int]](
        IndexedSeq[Int](0,13,21,22),
        IndexedSeq[Int](13,0,12,13),
        IndexedSeq[Int](21,12,0,13),
        IndexedSeq[Int](22,13,13,0)
    ))
    val d = getMutableDistanceMatrix(D)

    val a = additivePhylogeny(d,n)

    a.weight((0,4)) should be (11)
  }
}
