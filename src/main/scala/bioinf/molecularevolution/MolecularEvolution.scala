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

import bioinf.Input.finalizeMatrix
import bioinf.mutations.Mutations.AdjacencyList

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.math.max

object MolecularEvolution {

  /**
    * let dist be a |V| × |V| array of minimum distances initialized to ∞ (infinity)
    * let next be a |V| × |V| array of vertex indices initialized to null

    * procedure FloydWarshallWithPathReconstruction()
    *   for each edge (u,v)
    *     dist[u][v] ← w(u,v)  // the weight of the edge (u,v)
    *     next[u][v] ← v
    *   for k from 1 to |V| // standard Floyd-Warshall implementation
    *     for i from 1 to |V|
    *       for j from 1 to |V|
    *         if dist[i][k] + dist[k][j] < dist[i][j] then
    *           dist[i][j] ← dist[i][k] + dist[k][j]
    *           next[i][j] ← next[i][k]
    */
  def floydWarshallAlgorithm(x: AdjacencyList): FloydWarshall = {
    val n = max(x.v.max,x.w.max) + 1
    val dist = Array.fill[Int](n,n){x.w.sum}
    val next = Array.fill[Int](n,n){-1}

    for (i <- x.v.indices){
      val u = x.v(i)
      val v = x.w(i)
      val w = x.label(i)
      dist(u)(v) = w // the weight of the edge (u,v)
      next(u)(v) = v // next node
    }

    for (k <- dist.indices) {
      for (i <- dist.indices) { // standard Floyd-Warshall implementation
        for (j <- dist.indices) {
          if (dist(i)(k) + dist(k)(j) < dist(i)(j)) {
            dist(i)(j) = dist(i)(k) + dist(k)(j)
            next(i)(j) = next(i)(k)
          }
        }
      }
    }

    FloydWarshall(dist.map(_.toIndexedSeq).toIndexedSeq,next.map(_.toIndexedSeq).toIndexedSeq)
  }
  case class FloydWarshall(dist: IndexedSeq[IndexedSeq[Int]], next: IndexedSeq[IndexedSeq[Int]]) {
    def n = dist.length
    def apply(i: Int, j: Int) = dist(i)(j) // get distance between nodes
  }

  /**
    * procedure Path(u, v)
    *   if next[u][v] = null then
    *     return []
    *   path = [u]
    *     while u ≠ v
    *     u ← next[u][v]
    *     path.append(u)
    *   return path
    */
  def floydWarshallPathReconstruction(u: Int, v: Int, fw: FloydWarshall): IndexedSeq[Int] = {
    require (u < fw.n, u.toString + " was not less than " + fw.n.toString)
    require (v < fw.n, v.toString + " was not less than " + fw.n.toString)
    val path = new mutable.ArrayBuffer[Int]()
    if (fw.next(u)(v) >= 0) {
      path += u
      var currentNode = u
      while (currentNode != v) {
        currentNode = fw.next(currentNode)(v)
        path += currentNode
      }
    }
    path.result().toIndexedSeq
  }

  /**
    * In this chapter, we define the length of a path in a tree as the sum of the lengths of its edges (rather than the number of edges on the path). As a result, the evolutionary distance between two present-day species corresponding to leaves i and j in a tree T is equal to the length of the unique path connecting i and j, denoted di,j(T).

    * Distances Between Leaves Problem: Compute the distances between leaves in a weighted tree.
    * Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
    * Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.

    * CODE CHALLENGE: Solve the Distances Between Leaves Problem. The tree is given as an adjacency list of a graph whose leaves are integers between 0 and n - 1; the notation a->b:c means that node a is connected to node b by an edge of weight c. The matrix you return should be space-separated.

    * Sample Input:
    * 4
    * 0->4:11
    * 1->4:2
    * 2->5:6
    * 3->5:7
    * 4->0:11
    * 4->1:2
    * 4->5:4
    * 5->4:4
    * 5->3:7
    * 5->2:6
    * Sample Output:
    * 0	13	21	22
    * 13	0	12	13
    * 21	12	0	13
    * 22	13	13	0
    */
  def distanceBetweenLeaves(n: Int, x: AdjacencyList): IndexedSeq[IndexedSeq[Int]] = {
    val fw = floydWarshallAlgorithm(x)
    val a0 = new Array[Array[Int]](n)
    for (i <- 0 until n) {
      val a1 = Array.fill[Int](n){0}
      for (j <- 0 until n) {
        a1(j) = getDistanceBetweenLeaves(i,j,x,fw)
      }
      a0(i) = a1
    }
    a0.map{_.toIndexedSeq}.toIndexedSeq
  }

  def getDistanceBetweenLeaves(startNode: Int, endNode: Int, x: AdjacencyList, fw: FloydWarshall): Int = {
    val path = floydWarshallPathReconstruction(startNode,endNode,fw)
    var distance = 0
    for (i <- 0 until path.length - 1){
      distance += x.weight((path(i),path(i + 1)))
    }
    distance
  }

  /**
    *  We now have an algorithm for solving the Limb Length Problem.
    *  For each j, we can compute LimbLength(j) by finding the minimum value of
    *  (Di,j + Dj,k - Di,k)/2 over all pairs of leaves i and k.

    *  CODE CHALLENGE: Solve the Limb Length Problem.
    *       Input: An integer n, followed by an integer j between 0 and n - 1, followed by
    *         a space-separated additive distance matrix D (whose elements are integers).
    *       Output: The limb length of the leaf in Tree(D) corresponding to row j
    *         of this distance matrix (use 0-based indexing).
    *
    *  Sample Input:
    *       4
    *       1
    *       0	13	21	22
    *       13	0	12	13
    *       21	12	0	13
    *       22	13	13	0
    *  Sample Output:
    *       2
    */
  def limbLength(n: Int, j: Int, d: DistanceMatrix): Int = {
    // n = # of leaves
    // j = id of leaf for which to solve length
    require(j >= 0)
    require(j < d.n, j.toString + " exceeds maximum node Id " + (d.n - 1).toString)
    var min = Int.MaxValue
    var i = 0
    var k = 0
    while (i < d.n) {
      while (k < d.n) {
        /**
        *       we are searching for leaves i and k which minimize length of edge m->j
        * i    /
        *  >m-----j
        * k
        * ((i->m + m->j) + (k->m + m->j)) - (i->k)
        *   = path1 + path2 + 2*(common path) - common path
        */
        val x = (d(i,j) + d(k,j) - d(i,k)) / 2
        if (i != j && k != j && x < min) min = x
        k += 1
      }
      i += 1
    }
    min
  }

  def limb(d: DistanceMatrix, n: Int): Int = limbLength(d.n, n, d)

  /**
    *  CODE CHALLENGE: Implement AdditivePhylogeny to solve the Distance-Based Phylogeny Problem.
    *  Input: An integer n followed by a space-separated n x n distance matrix.
    *  Output: A weighted adjacency list for the simple tree fitting this matrix.
    *
    *  Note on formatting: The adjacency list must have consecutive integer node labels starting from 0.
    *  The n leaves must be labeled 0, 1, ..., n - 1 in order of their appearance in the distance matrix.
    *  Labels for internal nodes may be labeled in any order but must start from n and increase consecutively.
    *
    *  Sample Input:
    *    4
    *    0 13 21 22
    *    13 0	12 13
    *    21 12 0 13
    *    22 13 13 0
    *  Sample Output:
    *    0->4:11
    *    1->4:2
    *    2->5:6
    *    3->5:7
    *    4->0:11
    *    4->1:2
    *    4->5:4
    *    5->4:4
    *    5->3:7
    *    5->2:6
    *
    *    An algorithm for distance-based phylogeny reconstruction
    *
    *    The preceding discussion results in the following recursive algorithm, called
    *    AdditivePhylogeny, for finding the simple tree fitting an
    *    n x n additive distance matrix D.
    *    We assume that you have already implemented a program Limb(D, j) that computes
    *    LimbLength(j) for a leaf j based on the distance matrix D.
    *    Rather than selecting an arbitrary leaf j from Tree(D) for trimming,
    *    AdditivePhylogeny selects leaf n (corresponding to the last row and column of D).
    *
    *    AdditivePhylogeny(D, n)
    *        if n = 2
    *            return the tree consisting of a single edge of length D1,2
    *        limbLength ← Limb(D, n)
    *        for j ← 1 to n - 1
    *            Dj,n ← Dj,n - limbLength
    *            Dn,j ← Dj,n
    *        (i,n,k) ← three leaves such that Di,k = Di,n + Dn,k
    *        x ← Di,n
    *        remove row n and column n from D
    *        T ← AdditivePhylogeny(D, n - 1)
    *        v ← the (potentially new) node in T at distance x from i on the path between i and k
    *        add leaf n back to T by creating a limb (v, n) of length limbLength
    *        return T
    */
  def additivePhylogeny(d: MutableDistanceMatrix, nLeaves: Int): AdjacencyList = {
    require(nLeaves > 0)
    require(nLeaves <= d.n)
    // leaf indexes [0,n-1]
    // internal nodes [n,m-1]
    var n = nLeaves
    if (nLeaves > 2) {
      while (n > 2) {
        val j = n - 1 // max leaf node Id, index of nth leaf in DistanceMatrix
        val len = d.limbLength(j, nLeaves) // get length of last leaf
        d.addLimbLength(len, nLeaves)
        val path = d.findNodesForParent(nLeaves) // three nodes such that Di,k = Di,n + Dn,k
        val i = path._1
        val k = path._2
        val x = d(i, j) // x ← Di,n
        d.remove(j) // remove the last leaf
        val T = additivePhylogeny(d, nLeaves - 1) // "simple Tree fitting DistanceMatrix"
        // add node between i,k at distance x and return Tree'
        addLeaf(i, k, x, j, T)
        n -= 1
      }
      val a = AdjacencyListBuffer()
      a += Edge(0,1,d(0,1))
      a += Edge(1,0,d(0,1))
      a.result()
    } else { // a 2x2 distance matrix contains only nodes 0 and 1
      val a = AdjacencyListBuffer()
      a += Edge(0,1,d(0,1))
      a += Edge(1,0,d(0,1))
      a.result()
    }
  }

  case class DistanceMatrix(x: IndexedSeq[IndexedSeq[Int]]) {
    require(x.length >= 2)
    require(x.head.length == n)
    def n = x.length
    def indices = x.indices
    def apply(i: Int, j: Int): Int = x(i)(j)
  }

  case class MutableDistanceMatrix(d: mutable.ArrayBuffer[mutable.ArrayBuffer[Int]]) {
    require(d.length >= 2)
    require(d.length == d.head.length)
    def n = d.length
    def indices = d.indices
    def apply(i: Int, j: Int): Int = d(i)(j)
    def remove(i: Int) = {
      d.remove(i) // remove row
      d.foreach(_.remove(i)) // remove col
    }
    def limbLength(j: Int, nLeaves: Int): Int = {
      require(nLeaves <= n)
      // j = id of leaf for which to solve length
      require(j >= 0)
      require(j < nLeaves, j.toString + " exceeds maximum node Id " + (nLeaves - 1).toString)
      var min = Int.MaxValue
      var i = 0
      var k = 0
      while (i < nLeaves) {
        while (k < nLeaves) {
          /**
            *       we are searching for leaves i and k which minimize length of edge m->j
            * i    /
            *  >m-----j
            * k
            * ((i->m + m->j) + (k->m + m->j)) - (i->k)
            *   = path1 + path2 + 2*(common path) - common path
            */
          if (i != j && k != j){ // only examine edges between distinct leaves
            val x = (d(i)(j) + d(k)(j) - d(i)(k)) / 2
            if (x < min) min = x
          }
          k += 1
        }
        i += 1
      }
      min
    }
    def addLimbLength(limbLen: Int, nLeaves: Int) = {
      val n = nLeaves - 1
      for (j <- 0 until nLeaves){
        val edgeLen = d(j)(n)
        val newLen = edgeLen - limbLen // subtract limb length from distance
        // apply to both instances of the edge within the distance matrix
        d(j)(n) = newLen
        d(n)(j) = newLen
      }
    }
    def findNodesForParent(parentNode: Int): (Int,Int) = {
      var i = 0
      var k = 0
      val j = parentNode
      while (i < j) {
        while (k < j) {
          /** find leaves i,k with path through j
            *   |  i->j + j->k  |
            *   i ----- j ----- k
            *   |     i->k      |
            */
          if (d(i)(k) == d(i)(j) + d(j)(k)) {
            return (i, k)
          } else {
            k += 1
          }
        }
        i += 1
      }
      assert(false)
      (i,k)
    }
    def result(): DistanceMatrix = DistanceMatrix(d.map(_.result().toIndexedSeq).toIndexedSeq)
  }



  /** Search for node j in Tree T with distance x from node i and lies on path i->k
    *
    */
  def findNodeAtDistanceOnPath(i: Int, k: Int, x: Int, T: AdjacencyList): Option[Int] = {
    require(i < T.maxNodeId)
    require(k < T.maxNodeId)
    val fw = floydWarshallAlgorithm(T)
    val path = floydWarshallPathReconstruction(i, k, fw)
    /** search for node j with distance d from node i (may not exist)
      * i ---- j ---- k
      */
    path.find{j => fw.dist(i)(j) == x && j != i && j != k}
  }

  def getMutableDistanceMatrix(D: DistanceMatrix): MutableDistanceMatrix = {
    val m = new mutable.ArrayBuffer[mutable.ArrayBuffer[Int]]()
    D.x.foreach{row =>
      val rowBuf = new ArrayBuffer[Int](initialSize = row.length)
      rowBuf ++= row
      m += rowBuf
    }
    MutableDistanceMatrix(m)
  }

  // Remove row L and column L from DistanceMatrix
  def trimDistanceMatrix(D: DistanceMatrix, L: Int): DistanceMatrix = {
    val d1 = Array.fill[Int](D.n - 1, D.n - 1){Int.MinValue}
    var i1 = 0
    var j1 = 0
    for (i <- d1.indices) {
      for (j <- d1.indices) {
        if (i >= L) i1 = i + 1
        if (j >= L) j1 = j + 1
        d1(i)(j) = D(i1,j1)
      }
    }
    DistanceMatrix(finalizeMatrix(d1))
  }

  case class Edge(v: Int, w: Int, l: Int) {
    require(v != w)
    require(l >= 0)
  }

  case class AdjacencyListBuffer(v: mutable.ArrayBuffer[Int] = new mutable.ArrayBuffer[Int](), w: mutable.ArrayBuffer[Int] = new mutable.ArrayBuffer[Int](), label: mutable.ArrayBuffer[Int] = new mutable.ArrayBuffer[Int]()) {
    def result(): AdjacencyList = AdjacencyList(v.result().toIndexedSeq, w.result().toIndexedSeq, label.result().toIndexedSeq)
    def indices = v.indices
    def +=(e: Edge): Unit = {
      -=(e)
      v += e.v
      w += e.w
      label += e.l
    }
    def -=(e: Edge): Unit = indices.filter(i => v(i) == e.v && w(i) == e.w).foreach(remove)
    private def remove(i: Int): Unit = {
      v.remove(i)
      w.remove(i)
      label.remove(i)
    }
  }

  def getAdjacencyListBuffer(T: AdjacencyList): AdjacencyListBuffer = {
    val v = mutable.ArrayBuffer[Int]() ++= T.v
    val w = mutable.ArrayBuffer[Int]() ++= T.w
    val l = mutable.ArrayBuffer[Int]() ++= T.label
    AdjacencyListBuffer(v,w,l)
  }

  def addLeaf(i: Int, k: Int, x: Int, L: Int, T: AdjacencyList): AdjacencyList = {
    // Find node at distance x on path between nodes i and k
    val fw = floydWarshallAlgorithm(T)
    val path = floydWarshallPathReconstruction(i, k, fw)

    /** search for node j with distance d from node i (may not exist)
      * i ----- j ----- k
      * |<--d-->|
      */
    val node = path.find{j => fw(i,j) == x && j != i && j != k}
    val buf = getAdjacencyListBuffer(T)

    if (node.isEmpty) {
      var r = x // initialize distance countdown variable
      var p = i

      path.filter(_ != i).foreach{q =>
        r -= fw(p, q) // subtract distance between current nodes from countdown

        /**   exact match, use p (node on path i->k)
          *   |<- x ->|
          *   i -...- p -...- k
          */
        if (r == 0) {
          buf += Edge(L, q, x)
          buf += Edge(q, L, x)
        } else if (r < 0) {
          /**   No node at required distance; insert new node
            *                  |<-r-->|
            *   |<----- x ---->|
            *   i -...- p ----------- q -...- k
            *           |<-x0->|<-x1->|
            *
            *   x1 + x2 = |p->q|
            */
          val j = if (T.maxNodeId < L) {L} else {T.maxNodeId + 1}
          val x0 = r + fw(p,q)
          val x1 = -r

          //Adding leaf with limbLength at node newNode between origin and destination

          /** Remove edge between the two path nodes
            * t=0   i -...- p --------- q -...- k
            * t=1   i -...- p           q -...- k
            */
          buf -= Edge(p, q, 1)
          buf -= Edge(q, p, 1)

          /** Add edges between new node and first path node
            * t=1   i -...- p --------- q -...- k
            * t=2   i -...- p --- j     q -...- k
            */
          buf += Edge(p, j, x0)
          buf += Edge(j, p, x0)

          /** Add edges between new node and second path node
            * t=2   i -...- p --------- q -...- k
            * t=3   i -...- p --- j --- q -...- k
            */
          buf += Edge(q, j, x1)
          buf += Edge(j, q, x1)

          /** Add edges between new node and new leaf
            * t=3   i -...- p --- j --- q -...- k
            *
            * t=4                    L
            *                       /
            *       i -...- p --- j --- q -...- k
            */
          buf += Edge(j, L, x)
          buf += Edge(L, j, x)
        } else { // still have distance remaining
          p = q // prepare to evaluate next node in path
        }
      }
    }
    buf.result()
  }

}
