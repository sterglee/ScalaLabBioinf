/**
  * Consider the following adjacency list of a DAG:

  * a -> b: 5
  * a -> c: 6
  * a -> d: 5
  * b -> c: 2
  * b -> f: 9
  * c -> e: 4
  * c -> f: 3
  * c -> g: 7
  * d -> e: 4
  * d -> f: 5
  * e -> g: 2
  * f -> g: 1

  * What is the longest path in this graph? Give your answer as a sequence of nodes separated by spaces. (Note: a, b, c, d, e, f, g is a topological order for this graph.)

  */

import bioinf.comparinggenomes.ComparingGenomes
import ComparingGenomes._
import bioinf.Input._
val s = """a -> b: 5
a -> c: 6
a -> d: 5
b -> c: 2
b -> f: 9
c -> e: 4
c -> f: 3
c -> g: 7
d -> e: 4
d -> f: 5
e -> g: 2
f -> g: 1""".replaceAll(" ","")
val dag = readDAG(removeSpaces(char2num(s)).split("\n").toIndexedSeq)
val a1 = longestPathInDAG(0,6,dag)
val answer = num2char(a1)
