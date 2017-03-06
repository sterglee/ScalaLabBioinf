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

import bioinf.mutations.Mutations
import Mutations.{SuffixTree, AdjacencyList}
import bioinf.Input._

object Answer {

  def printVals(x: IndexedSeq[Int]): String = {
    val sb = new StringBuilder(1024)
    x.foreach{s =>
      sb.append(s)
      sb.append(" ")
    }
    sb.mkString
  }

  def printSeq(x: IndexedSeq[String]): String = {
    val sb = new StringBuilder(1024)
    x.foreach{s =>
      sb.append(s)
      sb.append(" ")
    }
    sb.mkString
  }

  def printTuples(x: IndexedSeq[(String,String)]): String = {
    val sb = new StringBuilder(1024000)
    for (i <- 0 until x.length - 1) {
      val tuple = x(i)
      sb.append(tuple._1)
      sb.append(" -> ")
      sb.append(tuple._2)
      sb.append("\n")
    }
    val tuple = x(x.length - 1)
    sb.append(tuple._1)
    sb.append(" -> ")
    sb.append(tuple._2)
    sb.mkString
  }

  def printDeBrujinGraph(x: IndexedSeq[(String,IndexedSeq[String])]): String = {
    val sb = new StringBuilder(1024000)
    val n = x.size
    var i = 0
    var j = 0
    x.foreach{tuple =>
      sb.append(tuple._1)
      sb.append(" -> ")
      j = 0
      tuple._2.foreach{pattern =>
        if (j > 0) sb.append(",")
        sb.append(pattern)
        j += 1
      }
      i += 1
      if (i < n) sb.append("\n")
    }
    sb.mkString
  }

  def printMatrix(x: IndexedSeq[IndexedSeq[Int]], sep: String = "\t"): String = {
    val sb = new StringBuilder(x.length * x.head.length * 4)
    for (i <- x.indices) {
      for (j <- x(i).indices) {
        sb.append(x(i)(j))
        if (j != x(i).indices.last) sb.append(sep)
      }
      if (i != x.indices.last) sb.append("\n")
    }
    sb.mkString
  }

  def printAdjacencyList(x: AdjacencyList): String = {
    val sb = new StringBuilder(x.v.length * 12)
    val indices = x.v.indices.sortBy(i => x.v(i)) // print in sorted order
    for (i <- indices) {
      sb.append(x.v(i))
      sb.append("->")
      sb.append(x.w(i))
      sb.append(':')
      sb.append(int2label(x.label(i)))
      sb.append("\n")
    }
    sb.delete(sb.length - 1,sb.length)
    sb.result()
  }

  def printTreeAdjacencyList(x: AdjacencyList): String = {
    val sb = new StringBuilder(x.v.length * 12)
    val indices = x.v.indices.sortBy(i => x.v(i)) // print in sorted order
    for (i <- indices) {
      sb.append(x.v(i))
      sb.append("->")
      sb.append(x.w(i))
      sb.append(':')
      sb.append(x.label(i))
      sb.append("\n")
    }
    sb.delete(sb.length - 1,sb.length)
    sb.result()
  }

  def printSuffixTreeEdgeLabels(x: SuffixTree): String = {
      val sb = new StringBuilder(65536)
      for (i <- x.edges.v.indices) {
        sb.append(x.text.substring(x.edges.pos(i),x.edges.pos(i) + x.edges.len(i)))
        sb.append("\n")
      }
      sb.delete(sb.length - 1, sb.length).result()
  }
}
