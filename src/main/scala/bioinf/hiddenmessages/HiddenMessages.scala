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

package bioinf.hiddenmessages

import java.nio.CharBuffer

import scala.collection.mutable

object HiddenMessages {

  /**
    * CODE CHALLENGE: Implement PatternCount (reproduced below).
    * Input: Strings Text and Pattern.
    * Output: Count(Text, Pattern).

    * PatternCount(Text, Pattern)
    * count ← 0
    * for i ← 0 to |Text| − |Pattern|
    * if Text(i, |Pattern|) = Pattern
    * count ← count + 1
    * return count

    * Sample Input:
    * GCGCG
    * GCG

    * Sample Output:
    * 2
    */
  def patternCount(text:String, pattern:String): Int = {
    var count = 0
    for (i <- 0 to text.length - pattern.length) {
      if (text.substring(i,i + pattern.length) == pattern)
        count += 1
    }
    count
  }

  /**
    * CODE CHALLENGE: Solve the Frequent Words Problem.
    * Input: A string Text and an integer k.
    * Output: All most frequent k-mers in Text.

    * Sample Input:
    * ACGTTGCATGTCGCATGATGCATGAGAGCT
    * 4

    * Sample Output:
    * CATG GCAT
    */
  def frequentWords(text:String, k:Int): IndexedSeq[String] = {
    val kmerCounts = kmerFrequency(text,k)
    val maxCount = kmerCounts.values.max // n for most common kmer(s)
    kmerCounts.filter(_._2 == maxCount).keys.toIndexedSeq.sorted // collection of most common kmers
  }

  // extracted for use in other algorithms
  // This function examines text for substrings of length k
  // it returns a map which allows lookup of the number of times a given kmer was found in the text
  @inline
  private def kmerFrequency(text:String,k:Int): Map[String,Int] = {
    getKmers(text,k).groupBy(s => s).map(t => (t._1,t._2.length)) // Map of (kmer,n)
  }

  @inline
  def getKmers(text:String,k:Int): IndexedSeq[String] = {
    val kmers = new Array[String](text.length - k + 1) // preallocate for performance instead of using map
    for (i <- 0 to text.length - k) {
      kmers(i) = text.substring(i, i + k)
    }
    kmers.toIndexedSeq.distinct
  }

  @inline
  private def complement(p:Char): Char = {
    p.toUpper match {
      case 'A' => 'T'
      case 'C' => 'G'
      case 'G' => 'C'
      case 'T' => 'A'
      case _ => 'N'
    }
  }

  /**
    * Given a nucleotide p, we denote its complementary nucleotide as p. The reverse complement of a string Pattern = p1…pn is the string Pattern = pn … p1 formed by taking the complement of each nucleotide in Pattern, then reversing the resulting string. We will need the solution to the following problem throughout this chapter:

    * Reverse Complement Problem: Find the reverse complement of a DNA string.
    * Input: A DNA string Pattern.
    * Output: Pattern, the reverse complement of Pattern.

    * CODE CHALLENGE: Solve the Reverse Complement Problem.

    * Sample Input:
    * AAAACCCGGT

    * Sample Output:
    * ACCGGGTTTT

    */
  def reverseComplement(pattern: String): String = {
    val nucleotides = CharBuffer.allocate(pattern.length)
    for (i <- pattern.indices.reverse) {
      nucleotides.put(complement(pattern.charAt(i)))
    }
    new String(nucleotides.array())
  }

  /**
    * CODE CHALLENGE: Solve the Pattern Matching Problem.
    * Input: Two strings, Pattern and Genome.
    * Output: A collection of space-separated integers specifying all starting positions where Pattern appears
    * as a substring of Genome.

    * Sample Input:
    * ATAT
    * GATATATGCATATACTT

    * Sample Output:
    * 1 3 9
    */
  def patternMatch(pattern:String,genome:String): IndexedSeq[Int] = {
    (0 to genome.length - pattern.length).flatMap{i =>
      genome.substring(i,i + pattern.length) match {
        case s if s == pattern => Some(i)
        case _ => None
      }
    }
  }

  /**
    * CODE CHALLENGE: Solve the Clump Finding Problem (restated below).
    * You will need to make sure that your algorithm is efficient enough to handle a large dataset.

    * Clump Finding Problem: Find patterns forming clumps in a string.
    * Input: A string Genome, and integers k, L, and t.
    * k size of k-mer
    * L size of window
    * t frequency threshold
    * Output: All distinct k-mers forming (L, t)-clumps in Genome.

    * You can solve the Clump Finding Problem by applying your algorithm for the
    * Frequent Words Problem to each window of length L in Genome.

    * Sample Input:
    * CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
    * 5 50 4

    * Sample Output:
    * CGACA GAAGA
    */
  def clumpFindingNaive(genome:String,k:Int,L:Int,t:Int): IndexedSeq[String] = {
    // this implementation does not use dynamic programming and is expected to be slow for large input
    // 150 second runtime for |genome| = 9128 k = 9 L = 598 t = 19
    (0 to genome.length - L).flatMap{i => // combine results for each window
      kmerFrequency(genome,k) // get
        .filter(_._2 >= t).keys
    }.sorted.distinct
  }

  /**
    * Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
           *Input: A DNA string Genome.
           *Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).

      *CODE CHALLENGE: Solve the Minimum Skew Problem.

      *Sample Input:
           *TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT

      *Sample Output:
           *11 24

      *values of Skewi (GAGCCACCGCGATA) for i ranging from 0 to 14.

      *Sample Input:
           *CATGGGCATCGGCCATACGCC

      *Sample Output:
           *0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
    */
  def getSkew(genome: IndexedSeq[Char]): IndexedSeq[Int] = {
    val a = Array.fill[Int](genome.length + 1){0}
    for (i <- genome.indices) {
      if (genome(i) == 'C') {
        a(i+1) = a(i) - 1
      } else if (genome(i) == 'G') {
        a(i+1) = a(i) + 1
      } else {
        a(i+1) = a(i)
      }
    }
    a.toIndexedSeq
  }
  def minimumSkew(genome: IndexedSeq[Char]): IndexedSeq[Int] = {
    val skew = getSkew(genome)
    skew.indices.filter(skew(_) == skew.min)
  }

  /**
    * Hamming Distance Problem: Compute the Hamming distance between two strings.
    * Input: Two strings of equal length.
    * Output: The Hamming distance between these strings.

    * CODE CHALLENGE: Solve the Hamming Distance Problem.

    * Sample Input:
    * GGGCCGTTGGT
    * GGACCGTTGAC
    * Sample Output:
    * 3
    */
  def hammingDistance(string1: String, string2: String): Int = {
    val buf1 = CharBuffer.wrap(string1).asReadOnlyBuffer()
    val buf2 = CharBuffer.wrap(string2).asReadOnlyBuffer()
    var distance = 0
    while (buf1.hasRemaining) {
      if (buf1.get() != buf2.get()) distance += 1
    }
    distance
  }

  /**
    *  We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear with slight variations leads to the following generalization of the Pattern Matching Problem.

    * Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
    * Input: Strings Pattern and Text along with an integer d.
    * Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

    * CODE CHALLENGE: Solve the Approximate Pattern Matching Problem.

    * Sample Input:
    * ATTCTGGA
    * CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
    * 3

    * Sample Output:
    * 6 7 26 27
    */
  def approximateMatch(a: String, b: String, c: Int): IndexedSeq[Int] = {
    val ab = mutable.ArrayBuffer[Int]()

    for (i <- 0 to b.length - a.length) {
      if (hammingDistance(a,b.substring(i, i + a.length)) <= c) ab += i
    }
    ab.result().toIndexedSeq
  }

  /**
    * Computing Countd(Text, Pattern) simply requires us to compute the Hamming distance between Pattern and every
    * k-mer substring of Text, which is achieved by the following pseudocode.

    * ApproximatePatternCount(Text, Pattern, d)
    * count ← 0
    * for i ← 0 to |Text| − |Pattern|
    * Pattern′ ← Text(i , |Pattern|)
    * if HammingDistance(Pattern, Pattern′) ≤ d
    * count ← count + 1
    * return count

    * STOP and Think: What is the running time of ApproximatePatternCount?

    * CODE CHALLENGE: Implement ApproximatePatternCount.
    * Input: Strings Pattern and Text as well as an integer d.
    * Output: Countd(Text, Pattern).

    * Visit the code-graded problem!

    * Extra Dataset

    * Debug Datasets

    * Sample Input:
    * GAGG
    * TTTAGAGCCTTCAGAGG
    * 2
    * Sample Output:
    * 4
    */
  def approximatePatternCount(pattern:String, text: String, d: Int): Int = {
    var count = 0
    for (i <- 0 to text.length - pattern.length) {
      val p1 = text.substring(i, i + pattern.length)
      if (hammingDistance(pattern, p1) <= d) {
        count += 1
      }
    }
    count
  }

  /**
    * A most frequent k-mer with up to d mismatches in Text is simply a string Pattern maximizing
    * Countd(Text, Pattern) among all k-mers. Note that Pattern does not need to actually appear as a
    * substring of Text; for example, as we already saw, AAAAA is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG, even though it does not appear exactly in this string. Keep this in mind while solving the following problem.

    * Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
    * Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
    * Output: All most frequent k-mers with up to d mismatches in Text.

    * One way to solve the above problem is to generate all 4k k-mers Pattern, compute
    * ApproximatePatternCount(Text, Pattern, d) for each k-mer Pattern, and then output k-mers with the
    * maximum number of approximate occurrences. This is an inefficient approach in practice,
    * since many of the 4k k-mers that this method analyzes should not be considered because neither
    * they nor their mutated versions (with up to d mismatches) appear in Text. Check out Charging Station:
    * Solving the Frequent Words with Mismatches Problem to learn about a better approach that avoids analyzing such hopeless k-mers.

    * CODE CHALLENGE: Solve the Frequent Words with Mismatches Problem.

    * Sample Input:
    * ACGTTGCATGTCGCATGATGCATGAGAGCT
    * 4 1

    * Sample Output:
    * GATG ATGC ATGT
    */
  def frequentWordsMismatch(text:String, k: Int, d: Int): IndexedSeq[String] = {
    val kmers = allKmers(k)
    val counts = Array.fill[Int](kmers.length){0}
    kmers.indices.foreach{i => counts(i) = approximatePatternCount(kmers(i),text,d)}
    val max = counts.max
    kmers.zip(counts).filter(_._2 == max).unzip._1
  }

  /**
    * We now redefine the Frequent Words Problem to account for both mismatches and reverse complements. Recall that Pattern refers to the reverse complement of Pattern.

      Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a string.
            Input: A DNA string Text as well as integers k and d.
            Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Pattern)
            over all possible k-mers.

      CODE CHALLENGE: Solve the Frequent Words with Mismatches and Reverse Complements Problem.

      Sample Input:
           ACGTTGCATGTCGCATGATGCATGAGAGCT
           4 1

      Sample Output:
           ATGT ACAT
    */
  def frequentWordsMismatchWithReverseComplement(text:String, k: Int, d: Int): IndexedSeq[String] = {
    val kmers = allKmers(k)
    val counts = Array.fill[Int](kmers.length){0}
    kmers.indices.foreach{i =>
      counts(i) = approximatePatternCount(kmers(i),text,d) + approximatePatternCount(reverseComplement(kmers(i)),text,d)
    }
    val max = counts.max
    kmers.zip(counts).filter{x => x._2 == max}.unzip._1
  }

  val NUCLEOTIDES = IndexedSeq[Char]('A','C','G','T')

  /**
    * Neighbors(Pattern, d)
    * if d = 0
    * return {Pattern}
    * if |Pattern| = 1
    * return {A, C, G, T}
    * Neighborhood ← an empty set
    * SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
    * for each string Text from SuffixNeighbors
    * if HammingDistance(Suffix(Pattern), Text) < d
    * for each nucleotide x
    * add x • Text to Neighborhood
    * else
    * add FirstSymbol(Pattern) • Text to Neighborhood
    * return Neighborhood
    *
    * */
  def neighbors(pattern: String, d: Int): IndexedSeq[String] = {
    val neighborhood = mutable.ArrayBuffer[String]()


    if (d == 0) {
      neighborhood += pattern
    } else if (pattern.length == 1) {
      NUCLEOTIDES.foreach{c => neighborhood += c.toString}
    } else {
      val chars = pattern.toCharArray
      val gen = chars.length
      var i = 0
      var j = 0
      var k = 0
      while (k < chars.length) {
        chars(i) = NUCLEOTIDES(j)
        neighborhood += new String(chars)
        if (j >= NUCLEOTIDES.length) {
          j = 0
          i += 1
        }
        if (i > k) {
          k += 1
        }
      }
    }
    neighborhood.result().toIndexedSeq
  }


  def pow(k: Int): Int = {
    require(k >= 0)
    var x = 1
    var _k = k
    while (_k > 0) {
      x *= 4
      _k -= 1
    }
    x
  }

  def allKmers(k: Int): IndexedSeq[String] = {
    val n = pow(k)
    val neighborhood = new Array[String](n)
    val chars = Array.fill[Char](k){'A'}
    val seq = Array.fill[Int](k){0}
    var m = 0

    def mutate(x: Array[Int]): Unit = {
      var i = x.length - 1
      var continue = true
      while (i >= 0 && continue) {
        // Increment last nucleotide
        if (x(i) == 3) {
          x(i) = 0
          i -= 1 // carry to previous nucleotides
        } else {
          x(i) += 1
          continue = false // no need to carry
        }
      }
    }

    while (m < n) {
      chars.indices.foreach{i => chars(i) = NUCLEOTIDES(seq(i))}
      neighborhood(m) = new String(chars)
      m += 1
      mutate(seq)
    }
    neighborhood.toIndexedSeq
  }
}
