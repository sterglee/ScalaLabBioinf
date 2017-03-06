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

package bioinf.comparinggenomes

import bioinf.UnitSpec
import bioinf.comparinggenomes.ComparingGenomes._
import bioinf.comparinggenomes.ScoringFunctions._

class ComparingGenomesSpec extends UnitSpec{
  "ComparingGenomes" should "score global alignmeent correctly" in {
    scoreAlignment("AAA","AAA",fn_BLOSUM62) should be (12)
  }

  it should "calculate BLOSUM62 score " in {
    fn_BLOSUM62('Y','Y') should be (7)
  }

  it should "calculate global alignment" in {
    val v = "PLEASANTLY"
    val w = "MEANLY"
    val sigma = 5
    val res = globalAlignment(v,w,sigma,fn_BLOSUM62)
    scoreGlobalAlignment(res) should be (8)
  }

  it should "calculate PAM250 score " in {
    fn_PAM250('Y','Y') should be (10)
  }

  it should "calculate local alignment" in {
    val v = "MEANLY"
    val w = "PENALTY"
    val sigma = 5
    val res = localAlignment(v,w,sigma,fn_PAM250)
    scoreLocalAlignment(res) should be (15)
  }

  it should "calculate affine gap alignment" in {
    val v = "PRTEINS"
    val w = "PRTWPSEIN"
    val sigma = 11
    val eta = 1
    val res = alignmentWithAffineGap(v, w, sigma, eta, fn_BLOSUM62)
    val score = scoreAffineGapAlignment(res)
    res.score should be (score)
    score should be (8)

  }

  it should "calculate longer affine gap alignment" in {
    val v = "KRYINALREEAYHCNNIHLFARCDDQRDNNYTQCTGYMGGVYYKWQFLIIQLYLCHSKVYAMSQMVVTPLRVTMYIV"
    val w = "KRALRNNIHLFARCDDPRDNNYTACTGYMGDVYYKWQFMIIHLYLCHSFQVYAMSQMVVEPLRVTMYEV"
    val sigma = 11
    val eta = 1
    val res = alignmentWithAffineGap(v, w, sigma, eta, fn_BLOSUM62)
    val score = scoreAffineGapAlignment(res)
    res.score should be (score)
    score should be (288)
    val print = res.print.split(System.lineSeparator())
    print(1) should be ("KRYINALREEAYHCNNIHLFARCDDQRDNNYTQCTGYMGGVYYKWQFLIIQLYLCHS-KVYAMSQMVVTPLRVTMYIV")
    print(2) should be ("KR---ALR------NNIHLFARCDDPRDNNYTACTGYMGDVYYKWQFMIIHLYLCHSFQVYAMSQMVVEPLRVTMYEV")
  }

  it should "align multiple sequences" in {
    val v = "AACGGGCAG"
    val w = "GGAAACAGAC"
    val u = "TAGCGGAAG"
    val res = multipleLongestCommonSubsequence(v,w,u,fn_AllEqual)
    val x = multipleLCS_Backtrack(res)
    val a = x.split(System.lineSeparator())
    a(0) should be ("4")
    a(1) should be ("A-----A-CGG-GC-AG---")
    a(2) should be ("-GGAA-A-C--AG--A-C--")
    a(3) should be ("-----TAGC---G-GA--AG")
  }

  it should "correctly implement greedy sorting of permutations" in {
    val p = IndexedSeq[Int](-3,4,1,5,-2)

    val buf = new Array[Int](p.length)
    reverse(0,0,p,buf).toIndexedSeq should be (IndexedSeq[Int](3,4,1,5,-2))
    reverse(0,1,p,buf).toIndexedSeq should be (IndexedSeq[Int](-4,3,1,5,-2))
    reverse(0,4,p,buf).toIndexedSeq should be (IndexedSeq[Int](2,-5,-1,-4,3))

    val res = greedySorting(p)
    val it = res.seq.iterator
    it.next() should be (IndexedSeq[Int](-1,-4,3,5,-2))
    it.next() should be (IndexedSeq[Int](1,-4,3,5,-2))
    it.next() should be (IndexedSeq[Int](1,2,-5,-3,4))
    it.next() should be (IndexedSeq[Int](1,2,3,5,4))
    it.next() should be (IndexedSeq[Int](1,2,3,-4,-5))
    it.next() should be (IndexedSeq[Int](1,2,3,4,-5))
    it.next() should be (IndexedSeq[Int](1,2,3,4,5))
    res.seq.length should be (7)
  }
}
