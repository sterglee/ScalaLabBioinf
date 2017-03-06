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

package bioinf.assembly

import bioinf.hiddenmessages.HiddenMessages

object Assembly {
  /**
    * CODE CHALLENGE: Solve the String Composition Problem.
    *Input: An integer k and a string Text.
      *Output: Compositionk(Text) (the k-mers can be provided in any order).

    *Sample Input:
      *5
      *CAATCCAAC
    *Sample Output:
      *CAATC
      *AATCC
      *ATCCA
      *TCCAA
      *CCAAC
  */
  def composition(k:Int,text:String): IndexedSeq[String] = {
    val kmers = new Array[String](text.length - k + 1) // preallocate for performance
    for (i <- 0 to text.length - k) {
      kmers(i) = text.substring(i,i + k)
    }
    kmers.distinct.sorted
  }

  /**
    * CODE CHALLENGE: Solve the String Spelled by a Genome Path Problem.

    * Sample Input:
    * ACCGA
    * CCGAA
    * CGAAG
    * GAAGC
    * AAGCT
    * Sample Output:
    * ACCGAAGCT
    */
  def genomePath(pathElems:IndexedSeq[String]): String = {
    assert(pathElems.length > 1)

    val seq = new Array[Char](pathElems.head.length + pathElems.length - 1)
    val it = pathElems.iterator
    val first = it.next()
    assert(first.length > 0)

    for (i <- seq.indices) {
      if (i < first.length) {
        seq(i) = first.charAt(i)
      } else {
        seq(i) = it.next().last
      }
    }
    new String(seq)
  }

  /**
    * CODE CHALLENGE: Solve the Overlap Graph Problem (restated below).
    * Input: A collection Patterns of k-mers.
    * Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the edges in any order.)

    * Sample Input:
    * ATGCG
    * GCATG
    * CATGC
    * AGGCA
    * GGCAT
    * Sample Output:
    * CATGC -> ATGCG
    * GCATG -> CATGC
    * GGCAT -> GCATG
    * AGGCA -> GGCAT
    */
  def overlap(patterns:IndexedSeq[String]): IndexedSeq[(String,String)] = {
    val prefixMap = prefixes(patterns)

    patterns.sorted.distinct.flatMap{pattern =>
      val suffix = pattern.substring(1, pattern.length)
      prefixMap.get(suffix) match {
        case Some(edges) =>
          edges.map{i =>
            (pattern,patterns(i))
          }
        case _ =>
          IndexedSeq[(String,String)]()
      }
    }
  }

  @inline
  private def prefixes(patterns:IndexedSeq[String]): Map[String,IndexedSeq[Int]] = {
    patterns.indices.groupBy{i =>
      patterns(i).substring(0, patterns(i).length - 1)
    }
  }

  /**
    * CODE CHALLENGE: Solve the De Bruijn Graph from a String Problem.
    * Input: An integer k and a string Text.
    * Output: DeBruijnk(Text), in the form of an adjacency list.

    * Sample Input:
    * 4
    * AAGATTCTCTAAGA
    * Sample Output:
    * AAG -> AGA,AGA
    * AGA -> GAT
    * ATT -> TTC
    * CTA -> TAA
    * CTC -> TCT
    * GAT -> ATT
    * TAA -> AAG
    * TCT -> CTA,CTC
    * TTC -> TCT
    */
  def deBrujinGraph(k: Int,text: String): IndexedSeq[(String,IndexedSeq[String])] = {
    val kmers = HiddenMessages.getKmers(text,k)
    deBrujinGraphFromKmers(kmers)
  }

  /**
    * DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
    * Input: A collection of k-mers Patterns.
    * Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).

    * CODE CHALLENGE: Solve the de Bruijn Graph from k-mers Problem.

    * Sample Input:
    * GAGG
    * CAGG
    * GGGG
    * GGGA
    * CAGG
    * AGGG
    * GGAG
    * Sample Output:
    * AGG -> GGG
    * CAG -> AGG,AGG
    * GAG -> AGG
    * GGA -> GAG
    * GGG -> GGA,GGG
    */
  def deBrujinGraphFromKmers(kmers: IndexedSeq[String]): IndexedSeq[(String,IndexedSeq[String])] = {
    val k = kmers.head.length
    val index = kmers.indices.groupBy(i => kmers(i).substring(0,k - 1)).toIndexedSeq.sortBy(_._1) // locations of each prefix

    index.map{tuple =>
      val prefix = tuple._1
      val adjacency = tuple._2.map{i =>
        kmers(i).substring(1,k)
      }.sorted
      (prefix,adjacency)
    }
  }

  // assumes that kmers overlap by k-1 characters
  def linearString(kmers: IndexedSeq[String]): String = {

    val alignments = {
      // Map prefixes to indices of the input array
      val pfx = kmers.indices.map{i =>
        val kmer = kmers(i)
        val k = kmer.length
        (kmer.substring(0,k-1),i)
      }.toMap

      // Map suffixes to indices of the matching prefix
      kmers.map {kmer =>
        val k = kmer.length
        (kmer, pfx.getOrElse(kmer.substring(1,k), -1))
      }.toMap
    }

    // get index of first kmer
    val firstKmer = alignments.find(_._2 == -1).getOrElse((kmers(0),0))._1
    val k = firstKmer.length
    val sb = new StringBuilder(kmers.length + k)
    sb.append(firstKmer.substring(0,firstKmer.length - 1)) // store the prefix of the first kmer
    var j = 0
    var i = 0
    while (i < kmers.length*2 && j >= 0) {
      sb.append(kmers(j).substring(kmers(j).length - 1, kmers(j).length)) // store the last character of the first kmer
      j = alignments.getOrElse(kmers(j),-1)
      i += 1
    }
    sb.result()
  }
}
