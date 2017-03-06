/**
    Question 3
    Imagine a hypothetical world in which there are two amino acids, X and Z,
    having respective masses 2 and 3. How many linear peptides can be formed
    from these amino acids having mass equal to 25?
    (Remember that the order of amino acids matters.)

  */

def getPermutations(lst: IndexedSeq[Int], target: Int, with_replacement: Boolean = true) = {
  val x = if (with_replacement) 0 else 1
  import scala.collection.mutable.ArrayBuffer
  def permute(lst: IndexedSeq[Int], idx: Int, l: ArrayBuffer[Int], r: ArrayBuffer[IndexedSeq[Int]], t: Int): ArrayBuffer[IndexedSeq[Int]] = {
    l.sum match {
      case x if x == t =>
        r += l.result().toIndexedSeq // append the sequence to the buffer
        l.clear()
      case x if x > t => // invalid sequence
      case _ =>
        for (i <- idx until lst.length) {
          permute(lst, i + x, l ++ Array(lst(i)), r, t) // append the character
        }
    }
    r
  }
  permute(lst, 0, ArrayBuffer[Int](), ArrayBuffer[IndexedSeq[Int]](), target)
}

implicit class Combinations(n: Int) {
  private def fact(n: Int): Int = (1 to n).foldLeft(1)(_ * _)
  def ! = fact(n) // allows 10!
  def choose(k: Int): Int = fact(n) / (fact(n - k) * fact(k))
}
// n Permutations of a characters of type i
def perms(a: Int, b: Int): Int = {
  ((a + b)!) / ((a!) * (b!))
}
def perms(counts: IndexedSeq[Int]): Int = {
  ((counts.sum)!) / (counts.map{_!}.foldLeft(1)(_ * _))
}
val lst = IndexedSeq(2,3)
val target = 25
val combinations = getPermutations(lst,target).result().toIndexedSeq
val counts = combinations.map{c =>
      lst.indices.map{i =>
        c.filter(_ == lst(i)).length
      }
    }

//val a = IndexedSeq(11,8,5,2)
//val b = IndexedSeq(1,3,5,7)
//val n = a.indices.map{i => perms(a(i),b(i))}
val n = counts.indices.map{i => perms(counts(i))}
val total = n.sum
