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

import bioinf.comparinggenomes.ComparingGenomes._

object ScoringFunctions {

  def fn_BLOSUM62(v: Char, w: Char): Int = {
    val iv = aminoAcidId(v)
    val iw = aminoAcidId(w)
    BLOSUM62(iv)(iw)
  }

  def fn_PAM250(v: Char, w: Char): Int = {
    val iv = aminoAcidId(v)
    val iw = aminoAcidId(w)
    PAM250(iv)(iw)
  }

  def fn_AllEqual(c1: Char, c2: Char, c3: Char): Int = {
    if (c1 == c2 && c2 == c3) 1
    else 0
  }

  val aminoAcidId = {
    val chars = IndexedSeq[Char]('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
    chars.indices.map{i => (chars(i),i)}.toMap
  }

  val BLOSUM62 = IndexedSeq[IndexedSeq[Int]](
    IndexedSeq(4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2),
    IndexedSeq(0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2),
    IndexedSeq(-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3),
    IndexedSeq(-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2),
    IndexedSeq(-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3),
    IndexedSeq(0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3),
    IndexedSeq(-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2),
    IndexedSeq(-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1),
    IndexedSeq(-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2),
    IndexedSeq(-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1),
    IndexedSeq(-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1),
    IndexedSeq(-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2),
    IndexedSeq(-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3),
    IndexedSeq(-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1),
    IndexedSeq(-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2),
    IndexedSeq(1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2),
    IndexedSeq(0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2),
    IndexedSeq(0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1),
    IndexedSeq(-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2),
    IndexedSeq(-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7)
  )

  val PAM250 = IndexedSeq[IndexedSeq[Int]](
    IndexedSeq[Int](2,-2,0,0,-3,1,-1,-1,-1,-2,-1,0,1,0,-2,1,1,0,-6,-3),
    IndexedSeq[Int](-2,12,-5,-5,-4,-3,-3,-2,-5,-6,-5,-4,-3,-5,-4,0,-2,-2,-8,0),
    IndexedSeq[Int](0,-5,4,3,-6,1,1,-2,0,-4,-3,2,-1,2,-1,0,0,-2,-7,-4),
    IndexedSeq[Int](0,-5,3,4,-5,0,1,-2,0,-3,-2,1,-1,2,-1,0,0,-2,-7,-4),
    IndexedSeq[Int](-3,-4,-6,-5,9,-5,-2,1,-5,2,0,-3,-5,-5,-4,-3,-3,-1,0,7),
    IndexedSeq[Int](1,-3,1,0,-5,5,-2,-3,-2,-4,-3,0,0,-1,-3,1,0,-1,-7,-5),
    IndexedSeq[Int](-1,-3,1,1,-2,-2,6,-2,0,-2,-2,2,0,3,2,-1,-1,-2,-3,0),
    IndexedSeq[Int](-1,-2,-2,-2,1,-3,-2,5,-2,2,2,-2,-2,-2,-2,-1,0,4,-5,-1),
    IndexedSeq[Int](-1,-5,0,0,-5,-2,0,-2,5,-3,0,1,-1,1,3,0,0,-2,-3,-4),
    IndexedSeq[Int](-2,-6,-4,-3,2,-4,-2,2,-3,6,4,-3,-3,-2,-3,-3,-2,2,-2,-1),
    IndexedSeq[Int](-1,-5,-3,-2,0,-3,-2,2,0,4,6,-2,-2,-1,0,-2,-1,2,-4,-2),
    IndexedSeq[Int](0,-4,2,1,-3,0,2,-2,1,-3,-2,2,0,1,0,1,0,-2,-4,-2),
    IndexedSeq[Int](1,-3,-1,-1,-5,0,0,-2,-1,-3,-2,0,6,0,0,1,0,-1,-6,-5),
    IndexedSeq[Int](0,-5,2,2,-5,-1,3,-2,1,-2,-1,1,0,4,1,-1,-1,-2,-5,-4),
    IndexedSeq[Int](-2,-4,-1,-1,-4,-3,2,-2,3,-3,0,0,0,1,6,0,-1,-2,2,-4),
    IndexedSeq[Int](1,0,0,0,-3,1,-1,-1,0,-3,-2,1,1,-1,0,2,1,-1,-2,-3),
    IndexedSeq[Int](1,-2,0,0,-3,0,-1,0,0,-2,-1,0,0,-1,-1,1,3,0,-5,-3),
    IndexedSeq[Int](0,-2,-2,-2,-1,-1,-2,4,-2,2,2,-2,-1,-2,-2,-1,0,4,-6,-2),
    IndexedSeq[Int](-6,-8,-7,-7,0,-7,-3,-5,-3,-2,-4,-4,-6,-5,2,-2,-5,-6,17,0),
    IndexedSeq[Int](-3,0,-4,-4,7,-5,0,-1,-4,-1,-2,-2,-5,-4,-4,-3,-3,-2,0,10)
  )

  def scoreGlobalAlignment(x: GlobalAlignmentResult): Int = {
    val s = x.print.split(System.lineSeparator())
    val v = s(1)
    val w = s(2)
    scoreAlignment(v,w,x.scoringFunction)
  }

  def scoreAlignment(v:String, w:String, f: (Char,Char) => Int): Int = {
    var score = 0
    v.indices.foreach{i =>
      if (v.charAt(i) == '-' || w.charAt(i) == '-') score -= 5
      else {
        score += f(v.charAt(i), w.charAt(i))
      }
    }
    score
  }

  def scoreLocalAlignment(x: LocalAlignmentResult): Int = {
    val s = x.print.split(System.lineSeparator())
    val v = s(1)
    val w = s(2)
    scoreAlignment(v,w,x.scoringFunction)
  }

  def scoreAffineGapAlignment(x: AffineGapAlignmentResult): Int = {
    val s = x.print.split(System.lineSeparator())
    val v = s(1)
    val w = s(2)
    scoreAffineGapAlignment(v, w, x.sigma, x.eta, x.scoringFunction)
  }

  def scoreAffineGapAlignment(v:String, w:String, sigma: Int, eta: Int, f: (Char,Char) => Int): Int = {
    var score = 0
    var vGapLength = 0
    var wGapLength = 0

    v.indices.foreach{i =>
      if (v.charAt(i) == '-') {
        if (vGapLength == 0){score -= sigma}
        else {score -= eta}
        vGapLength += 1
        wGapLength = 0
      } else if (w.charAt(i) == '-') {
        if (wGapLength == 0){score -= sigma}
        else {score -= eta}
        wGapLength += 1
        vGapLength = 0
      } else {
        score += f(v.charAt(i), w.charAt(i))
        vGapLength = 0
        wGapLength = 0
      }
    }
    score
  }




}
