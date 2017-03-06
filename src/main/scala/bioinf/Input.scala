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

package bioinf

import bioinf.mutations.Mutations.AdjacencyList

object Input {
  def getInts(s:String): IndexedSeq[Int] = {
    s.split(" ").map(_.toInt).toIndexedSeq
  }

  def transpose(x: IndexedSeq[IndexedSeq[Int]]):IndexedSeq[IndexedSeq[Int]] = {
    x.head.indices.map{i =>
      x.indices.map{j =>
        x(j)(i)
      }
    }
  }

  // edge notation: node0->node1:weight
  def readDAG(x: IndexedSeq[String]): IndexedSeq[IndexedSeq[Int]] = {
    val edges = Array.fill[Int](x.length,3){0}
    for (i <- x.indices) {
      val s = x(i)
      val i0: Int = s.indexOf("->")
      val i1: Int = s.indexOf(":")
      edges(i)(0) = s.substring(0,i0).toInt
      edges(i)(1) = s.substring(i0 + 2,i1).toInt
      edges(i)(2) = s.substring(i1 + 1,s.length).toInt
    }
    edges.map{_.toIndexedSeq}.toIndexedSeq
  }

  // edge notation: node0->node1:weight
  def readAdjacencyList(x: IndexedSeq[String]): AdjacencyList = {
    val v = Array.fill[Int](x.length){0}
    val w = Array.fill[Int](x.length){0}
    val label = Array.fill[Int](x.length){0}

    for (i <- x.indices) {
      val s = x(i)
      val i0: Int = s.indexOf("->")
      val i1: Int = s.indexOf(":")
      v(i) = s.substring(0,i0).toInt
      w(i) = s.substring(i0 + 2,i1).toInt
      label(i) = s.substring(i1 + 1,s.length).toInt
    }
    AdjacencyList(v,w,label)
  }

  def readMatrix(x: IndexedSeq[String]): IndexedSeq[IndexedSeq[Int]] = {
    val n = x.length
    val m = x.head.split(" ").length

    val a = Array.fill[Int](n,m){-1}
    for (i <- 0 until n) {
      for (j <- 0 until m) {
        val row = x(i).split(" ")
        assert(row.length == m)
        a(i)(j) = row(j).toInt
      }
    }
    finalizeMatrix(a)
  }

  @inline
  def finalizeMatrix(x: Array[Array[Int]]): IndexedSeq[IndexedSeq[Int]] = x.map(_.toIndexedSeq).toIndexedSeq

  @inline
  def removeSpaces(s: String): String = {
    s.replaceAll(" ","")
  }

  @inline
  def label2int(c: Char): Int = {
    c match {
      case 'A' => 0
      case 'C' => 1
      case 'G' => 2
      case 'T' => 3
      case '$' => -1
      case _ => -2
    }
  }

  @inline
  def int2label(x: Int): Char = {
    x match {
      case 0 => 'A'
      case 1 => 'C'
      case 2 => 'G'
      case 3 => 'T'
      case -1 => '$'
      case _ => 'N'
    }
  }

  @inline
  def char2num(s: String): String = {
    s.toLowerCase
      .replace('a','0')
      .replace('b','1')
      .replace('c','2')
      .replace('d','3')
      .replace('e','4')
      .replace('f','5')
      .replace('g','6')
      .replace('h','7')
      .replace('i','8')
      .replace('j','9')
  }

  @inline
  def num2char(s: String): String = {
    s.replace('0','a')
      .replace('1','b')
      .replace('2','c')
      .replace('3','d')
      .replace('4','e')
      .replace('5','f')
      .replace('6','g')
      .replace('7','h')
      .replace('8','i')
      .replace('9','j')
  }

  /** Sample Input:
    * (-3 +4 +1 +5 -2)
    * Output:
    * IndexedSeq[Int](-3,4,1,5,-2)
    */
  def readPermutation(s: String): IndexedSeq[Int] = {
    s.substring(1, s.length - 1).split(" ").map{_.toInt}.toIndexedSeq
  }
}
