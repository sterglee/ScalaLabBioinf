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

import bioinf.comparinggenomes.ComparingGenomes
import ComparingGenomes._
import bioinf.Answer._
import scala.io.Source
val it = Source.fromFile("""C:\tmp\dataset_245_5.txt""").getLines()
val s = it.next()
val t = it.next()
val paths = printMatrix(LCSPaths(s,t))
val bt = printMatrix(LCSBacktrack(s,t))
val ag = printMatrix(alignmentGraph(s,t))
val answer = longestCommonSubsequence(s,t)
// write results to disk
import java.nio.file.{Paths, Files}
import java.nio.charset.StandardCharsets.UTF_8
Files.write(Paths.get("""c:\tmp\results.txt"""),answer.getBytes(UTF_8))
Files.write(Paths.get("""c:\tmp\results_backtrack.txt"""),bt.getBytes(UTF_8))
Files.write(Paths.get("""c:\tmp\results_paths.txt"""),paths.getBytes(UTF_8))
Files.write(Paths.get("""c:\tmp\results_alignment.txt"""),ag.getBytes(UTF_8))
