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

package bioinf.mutations

import java.nio.CharBuffer

import bioinf.Input._
import bioinf.molecularevolution.MolecularEvolution
import MolecularEvolution.Edge
import bioinf.hiddenmessages.HiddenMessages

import scala.math.max
import scala.collection.mutable

object Mutations {
  /**
    * CODE CHALLENGE: Solve the Trie Construction Problem.
    * Input: A collection of strings Patterns.
    * Output: The adjacency list corresponding to Trie(Patterns), in the following format.
    * If Trie(Patterns) has n nodes, first label the root with 0 and then label the remaining nodes with
    * the integers 1 through n - 1 in any order you like. Each edge of the adjacency list of
    * Trie(Patterns) will be encoded by a triple: the first two members of the triple must be the
    * integers labeling the initial and terminal nodes of the edge, respectively; the third member
    * of the triple must be the symbol labeling the edge.

    * Sample Input:
    * ATAGA
    * ATC
    * GAT
    * Sample Output:
    * 0->1:A
    * 1->2:T
    * 2->3:A
    * 3->4:G
    * 4->5:A
    * 2->6:C
    * 0->7:G
    * 7->8:A
    * 8->9:T
    */
  def trieConstruction(patterns: IndexedSeq[String]): Trie = {
    val len = patterns.length * patterns.map {
      _.length
    }.max

    val node0 = Array.fill[Int](len) {
      -1
    }
    val node1 = Array.fill[Int](len) {
      -1
    }
    val label = Array.fill[Int](len) {
      -1
    }
    val edges = mutable.Map[(Int, Int), Int]() // edge pointers for (nodeId,label) tuples
    var i = 0 // edge pointer
    var j = 0 // node pointer

    patterns.foreach { pattern =>
      var currentNode = 0 // begin at ROOT node
      pattern.map {
        label2int
      }.foreach { currentSymbol =>
        val edge = (currentNode, currentSymbol)
        edges.get(edge) match {
          case Some(idx) => // edge exists
            currentNode = node1(idx)
          case _ => // new edge
            j += 1 // increment node pointer
            node0(i) = currentNode // write initial node to adjacency list
            node1(i) = j // write terminal node to adjacency list
            label(i) = currentSymbol // write label to adjacency list
            currentNode = j // continue from node j for next symbol
            edges += ((edge, i)) // add edge to index
            i += 1 // increment edge pointer
        }
      }

    }
    val data = IndexedSeq(node0, node1, label).map(_.slice(0, i).toIndexedSeq)
    Trie(AdjacencyList(data(0), data(1), data(2)))
  }

  case class AdjacencyList(v: IndexedSeq[Int], w: IndexedSeq[Int], label: IndexedSeq[Int]) {
    lazy val maxNodeId = max(v.max,w.max)
    lazy val out = v.indices.groupBy(v)
    lazy val weight = v.indices.map{i => ((v(i),w(i)),label(i))}.toMap
  }

  /**
    * CODE CHALLENGE: Implement TRIEMATCHING to solve the Multiple Pattern Matching Problem.
    * Input: A string Text and a collection of strings Patterns.
    * Output: All starting positions in Text where a string from Patterns appears as a substring.

    * TRIEMATCHING(Text, Trie)
    * while Text is nonempty
    * PREFIXTRIEMATCHING(Text, Trie)
    * remove first symbol from Text

    * PREFIXTRIEMATCHING(Text, Trie)
    * symbol ← first letter of Text
    * v ← root of Trie
    * while forever
    * if v is a leaf in Trie
    * return the pattern spelled by the path from the root to v
    * else if there is an edge (v, w) in Trie labeled by symbol
    * symbol ← next letter of Text
    * v ← w
    * else
    * output "no matches found"
    * return

    * Sample Input:
    * AATCGGGTTCAATCGGGGT
    * ATCG
    * GGGT
    * Sample Output:
    * 1 4 11 15
    */
  def trieMatching(text: String, patterns: IndexedSeq[String]): IndexedSeq[Int] = {
    val buf = CharBuffer.wrap(text).asReadOnlyBuffer()
    val trie = trieConstruction2(patterns)
    val matchIndices = mutable.ArrayBuffer[Int]()

    var symbol = 0

    while (buf.hasRemaining) {
      symbol = label2int(buf.get())
      buf.mark()
      var v = 0
      var continue = true
      while (continue) {
        val edge = trie.edges0.get((v, symbol))
        val leaf = trie.edges0.get((v, -1)) // pattern termination
        if (leaf.nonEmpty) {
          // match is complete
          buf.reset()
          matchIndices += buf.position() - 1
          continue = false
        } else if (edge.nonEmpty && buf.hasRemaining) {
          symbol = label2int(buf.get()) // next letter of text
          v = trie.adjacencyList.w(edge.get) // next node in trie
        } else {
          // no matches found
          buf.reset()
          continue = false
        }
      }
    }
    matchIndices.result().toIndexedSeq
  }

  case class Trie(adjacencyList: AdjacencyList) {
    lazy val leaves = adjacencyList.w.toSet.diff(adjacencyList.v.toSet)
    lazy val edges = adjacencyList.label.indices.map { i => ((adjacencyList.v(i), adjacencyList.w(i)), i) }.toMap
    lazy val edges0 = adjacencyList.label.indices.map { i => ((adjacencyList.v(i), adjacencyList.label(i)), i) }.toMap
    lazy val edges1 = adjacencyList.label.indices.map { i => ((adjacencyList.w(i), adjacencyList.label(i)), i) }.toMap
  }

  /** This trie construction method appends '$' to each pattern to indicate
    * This allows successful resolution of patterns which are nested within each other
    * Example in the book is 'pan' and 'pantry'
    * With terminator characters, they become 'pan$' and 'pantry$' such that 'pan' can be found in the trie
    */
  def trieConstruction2(patterns: IndexedSeq[String]): Trie = {
    trieConstruction(patterns.map {
      _ + "$"
    })
  }

  /**
    * CODE CHALLENGE: Solve the Suffix Tree Construction Problem.
    * Input: A string Text.
    * Output: The edge labels of SuffixTree(Text). You may return these strings in any order.

    * Sample Input:
    * ATAAATG$
    * Sample Output:
    * AAATG$
    * G$
    * T
    * ATG$
    * TG$
    * A
    * A
    * AAATG$
    * G$
    * T
    * G$
    * $

    * MODIFIEDSUFFIXTRIECONSTRUCTION(Text)
    * Trie ← a graph consisting of a single node root
    * for i ← 0 to |Text| - 1
    * currentNode ← root
    * for j ← i to |Text| - 1
    * currentSymbol ← j-th symbol of Text
    * if there is an outgoing edge from currentNode labeled by currentSymbol
    * currentNode ← ending node of this edge
    * else
    * add a new node newNode to Trie
    * add an edge newEdge connecting currentNode to newNode in Trie
    * Symbol(newEdge) ← currentSymbol
    * Position(newEdge) ← j
    * currentNode ← newNode
    * if currentNode is a leaf in Trie
    * assign label i to this leaf
    * return Trie

    *
    */
  def suffixTrieConstruction(text: String): SuffixTrie = {
    val v = mutable.ArrayBuffer[Int]()
    val w = mutable.ArrayBuffer[Int]()
    val pos = mutable.ArrayBuffer[Int]()
    val label = mutable.ArrayBuffer[Int]()
    var edgeId = 0 // edge pointer
    var nodeId = 0 // node pointer
    val edges = mutable.Map[(Int, Int), Int]() // edge pointers for (nodeId,label) tuples
    val positionLabel = mutable.Map[Int, Int]() // index pointers for edgeId -> position
    var currentNode = 0
    var symbol = 0
    var edge = (0, 0)
    val TERMINATOR_SYMBOL = -1

    for (i <- 0 until text.length) {
      currentNode = 0
      for (j <- i until text.length) {
        symbol = label2int(text.charAt(j))
        edge = (currentNode, symbol)
        edges.get(edge) match {
          case Some(idx) =>
            currentNode = w(idx)
          case _ =>
            // add a new node
            v += currentNode
            nodeId += 1
            w += nodeId

            // add an edge
            edges += ((edge, edgeId))
            pos += j
            label += symbol
            edgeId += 1

            currentNode = nodeId
        }

      }
      if (symbol == TERMINATOR_SYMBOL) {
        positionLabel += ((edges.getOrElse(edge, 0), i))
      }
    }
    val data = IndexedSeq(v, w, label, pos).map {
      _.result().toIndexedSeq
    }
    SuffixTrie(AdjacencyList(data(0), data(1), data(2)), data(3), positionLabel.toMap)
  }

  case class SuffixTrie(adjacencyList: AdjacencyList, position: IndexedSeq[Int], label: Map[Int, Int]) {
    def indices = adjacencyList.v.indices

    def length = adjacencyList.v.length

    lazy val outEdges = indices.groupBy(i => adjacencyList.v(i))
  }

  case class Path(v: Int, w: Int, pos: Int, len: Int)

  /**
    * MODIFIEDSUFFIXTREECONSTRUCTION(Text)
    * Trie ← MODIFIEDSUFFIXTRIECONSTRUCTION
    * for each non-branching path Path in Trie
    * substitute Path by a single edge e connecting the first and last nodes of Path
    * Position(e) ← Position(first edge of Path)
    * Length(e) ← number of edges of Path
    * return Trie
    */
  def suffixTreeConstruction(text: String): SuffixTree = {
    val text1 = if (text.last == '$') text else text + "$"
    val t = suffixTrieConstruction(text1)

    val stack = mutable.Stack[Int]() // used for DFS
    val v = mutable.ArrayBuffer[Int]()
    val w = mutable.ArrayBuffer[Int]()
    val pos = mutable.ArrayBuffer[Int]()
    val len = mutable.ArrayBuffer[Int]()

    val out = t.outEdges
    val currentNode = 0
    out.get(currentNode).foreach(_.foreach(stack.push))

    while (stack.nonEmpty) {
      val i = stack.pop()
      val startNode = t.adjacencyList.v(i)
      val position = t.position(i)

      var length = 1
      var currentNode = t.adjacencyList.w(i)
      var next: Option[Int] = t.outEdges.get(currentNode) match {
        case Some(ids) if ids.length == 1 => // Path
          Some(ids.head)
        case Some(ids) if ids.length > 1 => // Branch
          ids.foreach(stack.push)
          None
        case _ => // leaf
          None
      }
      while (next.nonEmpty) {
        currentNode = t.adjacencyList.w(next.get)
        length += 1
        next = t.outEdges.get(currentNode) match {
          case Some(ids) if ids.length == 1 => // Path
            Some(ids.head)
          case Some(ids) if ids.length > 1 => // Branch
            ids.foreach(stack.push)
            None
          case _ => // leaf
            None
        }
      }
      v += startNode
      w += currentNode
      pos += position
      len += length
    }

    SuffixTree(text1, Edges(v.toIndexedSeq, w.toIndexedSeq, pos.toIndexedSeq, len.toIndexedSeq))
  }

  case class Edges(v: IndexedSeq[Int], w: IndexedSeq[Int], pos: IndexedSeq[Int], len: IndexedSeq[Int])

  case class SuffixTree(text: String, edges: Edges) {
    lazy val out = edges.v.indices.groupBy(i => edges.v(i))
    lazy val in = edges.v.indices.map { i => (edges.w(i), i) }.toMap
  }

  // returns indices of edges terminating at deepest nodes
  def deepestEdges(s: SuffixTree): IndexedSeq[Int] = {
    val depth = Array.fill[Int](s.edges.v.length) {
      -1
    }
    val stack = mutable.Stack[Int]() // used for DFS

    var currentNode = 0
    var parentNode = -1
    var maxDepth = -1
    val deepEdgeIds = mutable.Set[Int]()

    s.out.get(currentNode).foreach(_.foreach(stack.push))
    while (stack.nonEmpty) {
      val i = stack.pop()
      parentNode = s.edges.v(i)
      currentNode = s.edges.w(i)
      val currentDepth = s.in.get(parentNode) match {
        case Some(idx) => depth(idx) + s.edges.len(i)
        case _ => s.edges.len(i)
      }
      depth(i) = currentDepth
      s.out.get(currentNode) match {
        case Some(indices) =>
          if (indices.length > 1) {
            // only return edges that branch
            if (currentDepth > maxDepth) {
              maxDepth = currentDepth
              deepEdgeIds.clear()
              deepEdgeIds += i
            } else if (currentDepth == maxDepth) {
              deepEdgeIds += i
            }
          }
          indices.foreach(stack.push)
        case _ =>
      }
    }
    deepEdgeIds.toIndexedSeq
  }

  // prints the string terminating at edge with given index
  def printSuffix(s: SuffixTree, startIdx: Int): String = {
    var parentNode = s.edges.v(startIdx)
    val sb = new StringBuilder(65536)
    val stack = mutable.Stack[Int]()
    stack.push(startIdx)
    while (stack.nonEmpty) {
      val i = stack.pop()
      parentNode = s.edges.v(i)
      s.text
        .substring(s.edges.pos(i), s.edges.pos(i) + s.edges.len(i))
        .reverseIterator
        .foreach(sb.append)
      s.in.get(parentNode).foreach(stack.push)
    }
    sb.reverseContents().mkString
  }

  val ROOT = 0
  val GREY = -1
  val RED = 1
  val BLUE = 2
  val PURPLE = 3

  // returns colors of all edges
  def getNodeColors(s: SuffixTree): IndexedSeq[Int] = {
    val color = Array.fill[Int](s.edges.v.length) {
      GREY
    }
    val stack = mutable.Stack[Int]() // used for DFS
    val stack2 = mutable.Stack[Int]() // used for bottom-up traversal
    val split = s.text.indexOf("#")

    // Depth-first search
    s.out.get(ROOT).foreach(_.foreach(stack.push))
    while (stack.nonEmpty) {
      val i = stack.pop()
      val currentNode = s.edges.w(i)
      s.out.get(currentNode) match {
        case Some(indices) =>
          indices.foreach(stack.push)
        case _ => // leaf
          if (s.edges.pos(i) > split - 1) {
            color(i) = RED
          } else {
            color(i) = BLUE
          }
          stack2.push(i)
      }
    }

    // Bottom-up traversal
    while (stack2.nonEmpty) {
      val i = stack2.pop()
      val c = color(i)

      // Only search children if color is Gray
      if (c == GREY) {
        val currentNode = s.edges.w(i)
        val children = s.out.get(currentNode)
        val newColor = children match {
          case Some(indices) =>
            val colors = indices.map(color)
            if (colors.contains(RED) && colors.contains(BLUE)) PURPLE
            else colors.head // Child color
          case _ =>
            Int.MinValue // should never happen
        }
        color(i) = newColor
      }

      val parentNode = s.edges.v(i)
      val parent = s.in.get(parentNode)
      parent match {
        case Some(idx) if color(idx) == GREY => stack2.push(idx)
        case _ => // root
      }
    }

    color.toIndexedSeq
  }

  // returns indices of edges terminating at deepest nodes
  def deepestColoredEdges(s: SuffixTree): IndexedSeq[Int] = {
    val depth = Array.fill[Int](s.edges.v.length) {
      -1
    }
    val stack = mutable.Stack[Int]() // used for DFS

    var currentNode = 0
    var parentNode = -1
    var maxDepth = -1
    val edges = mutable.Set[Int]()
    val color = getNodeColors(s)

    s.out.get(currentNode).foreach(_.foreach(stack.push))
    while (stack.nonEmpty) {
      val i = stack.pop()
      parentNode = s.edges.v(i)
      currentNode = s.edges.w(i)
      val currentDepth = s.in.get(parentNode) match {
        case Some(idx) => depth(idx) + s.edges.len(i)
        case _ => s.edges.len(i)
      }
      depth(i) = currentDepth
      s.out.get(currentNode) match {
        case Some(indices) =>
          if (indices.length > 1 && color(i) == PURPLE) {
            // only return edges with colored branch nodes
            if (currentDepth > maxDepth) {
              maxDepth = currentDepth
              edges.clear()
              edges += i
            } else if (currentDepth == maxDepth) {
              edges += i
            }
          }
          indices.foreach(stack.push)
        case _ =>
      }
    }
    edges.toIndexedSeq
  }

  /**
    * Shortest Non-Shared Substring Problem: Find the shortest substring of one string that does not appear in another string.
    * Input: Strings Text1 and Text2.
    * Output: The shortest substring of Text1 that does not appear in Text2.

    * CODE CHALLENGE: Solve the Shortest Non-Shared Substring Problem. (Multiple solutions may exist, in which case you may return any one.)

    * Sample Input:
    * CCAAGCTGCTAGAGG
    * CATGCTGGGCTGGCT

    * Sample Output:
    * AA
    */
  def shortestNonSharedSubstring(text1: String, text2: String): IndexedSeq[String] = {
    val s2 = suffixTreeConstruction(text2)
    val nonSharedSubstrings = mutable.ArrayBuffer[String]()
    var k = 1
    while (nonSharedSubstrings.isEmpty) {
      val kmers = HiddenMessages.getKmers(text1, k)
      kmers.foreach { s =>
        if (!patternExistsInSuffixTree(s, s2)) nonSharedSubstrings += s
      }
      k += 1
    }
    nonSharedSubstrings.result().sorted.toIndexedSeq
  }

  def patternExistsInSuffixTree(pattern: String, s: SuffixTree): Boolean = {
    val matchIndices = mutable.ArrayBuffer[Int]()
    val stack = mutable.Stack[Int]()
    var beginIndex = 0
    var endIndex = 1

    s.out.get(ROOT) match {
      case Some(ids) => ids.foreach { i =>
        if (s.text.charAt(s.edges.pos(i)) == pattern.charAt(0)) stack.push(i)
      }
      case _ =>
    }
    while (stack.nonEmpty) {
      val i = stack.pop()
      val pos = s.edges.pos(i)
      val compareLen = scala.math.min(s.edges.len(i), pattern.length - beginIndex)
      endIndex = beginIndex + compareLen
      val ps = pattern.substring(beginIndex, endIndex)
      val ns = s.text.substring(pos, pos + compareLen)
      beginIndex = endIndex
      if (ps == ns) {
        if (endIndex == pattern.length) matchIndices += pos
        else {
          s.out.get(s.edges.w(i)) match {
            case Some(ids) => ids.foreach { i =>
              if (s.text.charAt(s.edges.pos(i)) == pattern.charAt(beginIndex)) stack.push(i)
            }
            case _ =>
          }
        }
      }
    }
    if (matchIndices.nonEmpty) true
    else false
  }

  /**
    * Suffix Array Construction Problem: Construct the suffix array of a string.
    * Input: A string Text.
    * Output: SuffixArray(Text).

    * CODE CHALLENGE: Solve the Suffix Array Construction Problem.

    * Sample Input:
    * AACGATAGCGGTAGA$

    * Sample Output:
    * 15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5
    */
  def suffixArrayConstruction(text: String): IndexedSeq[Int] = {
    text.indices
      .map { i => (i, text.substring(i)) }
      .sortBy(_._2)
      .map {
        _._1
      }
  }

  /**
    * Burrows-Wheeler Transform Construction Problem: Construct the Burrows-Wheeler transform of a string.
    * Input: A string Text.
    * Output: BWT(Text).

    * CODE CHALLENGE: Solve the Burrows-Wheeler Transform Construction Problem.

    * Note: Although it is possible to construct the Burrows-Wheeler transform in O(|Text|) time and space, we do not expect you to implement such a fast algorithm. In other words, it is perfectly fine to produce BWT(Text) by first producing the complete Burrows-Wheeler matrix M(Text).

    * Sample Input:
    * GCGTGCCTGGTCA$
    * Sample Output:
    * ACTGGCT$TGCGGC*

    * Get all cyclic rotations
    * sort the rotations
    * take the last character of each

    */
  def burrowsWheelerTransformConstruction(text: String): String = {
    require(text.length > 1, "text length must be greater than 1")
    val cyclicRotations = new Array[String](text.length())
    cyclicRotations(0) = rotate(text)
    for (i <- 1 until cyclicRotations.length) {
      cyclicRotations(i) = rotate(cyclicRotations(i - 1))
    }
    val bwt = cyclicRotations.sorted.map {
      _.charAt(text.length - 1)
    }
    new String(bwt)
  }

  def rotate(chars: String): String = {
    require(chars.length > 1, "text length must be greater than 1")
    chars.substring(chars.length - 1, chars.length) + chars.substring(0, chars.length - 1)
  }

  /**
    * Inverse Burrows-Wheeler Transform Problem: Reconstruct a string from its Burrows-Wheeler transform.
    * Input: A string Transform (with a single "$" symbol).
    * Output: The string Text such that BWT(Text) = Transform.

    * CODE CHALLENGE: Solve the Inverse Burrows-Wheeler Transform Problem.

    * Sample Input:
    * TTCCTAACG$A
    * Sample Output:
    * TACATCACGT$
    */
  def inverseBurrowsWheelerTransform(transform: String): String = {
    val col1 = transform.toCharArray.sorted
    val col1Indices = getLastToFirstMapping(transform.toCharArray.toIndexedSeq)

    val buf = CharBuffer.allocate(col1.length) // used to output original sequence

    var i = transform.indexOf("$") // the first column at this index will contain first character of the original sequence
    var c = col1(i)
    buf.put(c) // output first character of original sequence
    while (c != '$') {
      // collect characters until terminator character is reached
      i = col1Indices(i) // locate the next character within first column
      c = col1(i)
      buf.put(c) // output a character
    }
    buf.flip()
    buf.toString
  }

  // Used to locate index of nth instance of each character in first column
  def getCol1IndicesByChar(bwt: IndexedSeq[Char]): IndexedSeq[Int] = {
    val col1 = bwt.sorted.toArray
    val col1Index = mutable.ArrayBuffer[Int]() // stores index within col1 for each character of BWT
    bwt.foreach { c =>
      val i = col1.indexOf(c)
      col1Index += i
      col1.update(i, 0) // prevent this index from being found again
    }
    col1Index.result().toIndexedSeq
  }

  /**
    * Given a symbol at position i of LastColumn, the Last-to-First mapping identifies this symbol’s
    * position in FirstColumn.
    */
  def getLastToFirstMapping(bwt: IndexedSeq[Char]): IndexedSeq[Int] = {
    val firstIndices = getCol1IndicesByChar(bwt)
    bwt.indices.map { i => firstIndices.indexOf(i) } // indices of first column corresponding to kth character of BWT
  }

  /**
    * CODE CHALLENGE: Implement BWMATCHING.
    * Input: A string BWT(Text), followed by a collection of Patterns.
    * Output: A list of integers, where the i-th integer corresponds to the number of substring
    * matches of the i-th member of Patterns in Text.

    * Extra Dataset

    * Sample Input:
    * TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC
    * CCT CAC GAG CAG ATC
    * Sample Output:
    * 2 1 1 0 1
    *
    * We are now ready to describe BWMATCHING, an algorithm that counts the total number of matches of Pattern in Text, where the only information that we are given is FirstColumn and LastColumn in addition to the Last-to-First mapping. The pointers top and bottom are updated by the green lines in the following pseudocode.

    * BWMATCHING(FirstColumn, LastColumn, Pattern, LastToFirst)
    * top ← 0
    * bottom ← |LastColumn| − 1
    * while top ≤ bottom
    * if Pattern is nonempty
    * symbol ← last letter in Pattern
    * remove last letter from Pattern
    * if positions from top to bottom in LastColumn contain an occurrence of symbol
    * topIndex ← first position of symbol among positions from top to bottom in LastColumn
    * bottomIndex ← last position of symbol among positions from top to bottom in LastColumn
    * top ← LastToFirst(topIndex)
    * bottom ← LastToFirst(bottomIndex)
    * else
    * return 0
    * else
    * return bottom − top + 1
    */
  //def bwMatching(bwt: String, patterns: IndexedSeq[String]): IndexedSeq[Int] = {IndexedSeq[Int]()}

  /**
    * Define FirstOccurrence(symbol) as the first position of symbol in FirstColumn.
    * If Text = "panamabananas$", then FirstColumn is "$aaaaaabmnnnps", and the array holding all values of
    * FirstOccurrence is [0, 1, 7, 8, 9, 11, 12], as shown below. For DNA strings of any length,
    * the array FirstOccurrence contains only five elements
    */
  def getFirstOccurrence(bwt: IndexedSeq[Char]): Map[Char,Int] = {
    val firstColumn = bwt.sorted
    bwt.distinct.sorted.map{c => (c,firstColumn.indexOf(c))}.toMap
  }

  /**
    * To improve BWMATCHING, we introduce a function Count(symbol, i, LastColumn), which returns the number of
    * occurrences of symbol in the first i positions of LastColumn.
    * For example, Count"n”(10, "smnpbnnaaaaa$a”) = 3, and Count"a”(4, "smnpbnnaaaaa$a”) = 0.
    * Below, we show arrays holding Countsymbol(i, "smnpbnnaaaaa$a") for every symbol occurring in "panamabananas$”.
    */
  def getCounts(symbol: Char, lastColumn: IndexedSeq[Char]): IndexedSeq[Int] = {
    val symbolCount = Array.fill[Int](lastColumn.length + 1) {0}
    for (i <- symbolCount.indices) {
      symbolCount(i) = getCount(symbol, i, lastColumn)
    }
    symbolCount.toIndexedSeq
  }

  @inline
  def getCount(symbol: Char, i: Int, lastColumn: IndexedSeq[Char]): Int = {
    var n = 0
    var j = 0
    while (j < i) {
      if (lastColumn(j) == symbol) n += 1
      j += 1
    }
    n
  }

  /**
    * BETTERBWMATCHING(FirstOccurrence, LastColumn, Pattern, Count)
    * top ← 0
    * bottom ← |LastColumn| − 1
    * while top ≤ bottom
    *   if Pattern is nonempty
    *     symbol ← last letter in Pattern
    *     remove last letter from Pattern
    *     if positions from top to bottom in LastColumn contain an occurrence of symbol
    *       top ← FirstOccurrence(symbol) + Countsymbol(top, LastColumn)
    *       bottom ← FirstOccurrence(symbol) + Countsymbol(bottom + 1, LastColumn) − 1
    *     else
    *       return 0
    *   else
    *     return bottom − top + 1
    */

  def betterBWMatching(firstOccurrence: Map[Char,Int], lastColumn: IndexedSeq[Char], pattern: IndexedSeq[Char], count: Map[Char,IndexedSeq[Int]]): Int = {
    var top = 0
    var bottom = lastColumn.length - 1
    val p = mutable.Stack[Char]()
    pattern.foreach{c => p.push(c)}
    var symbolOccurs = true
    while (p.nonEmpty && top <= bottom) {
      val symbol = p.pop()
      symbolOccurs = lastColumn.slice(top,bottom+1).contains(symbol)
      if (symbolOccurs) {
        top = firstOccurrence(symbol) + count(symbol)(top)
        bottom = firstOccurrence(symbol) + count(symbol)(bottom + 1) - 1
      } else return 0
    }
    bottom - top + 1
  }

  def betterBurrowsWheelerMatching(bwt: String, patterns: IndexedSeq[String]): IndexedSeq[Int] = {
    val bwtChars = bwt.toCharArray.toIndexedSeq
    val firstOccurrence = getFirstOccurrence(bwtChars)
    val count = bwtChars.distinct.map{c => (c, getCounts(c, bwtChars))}.toMap

    patterns.map{pattern =>
      betterBWMatching(firstOccurrence, bwtChars, pattern, count)
    }
  }
}
