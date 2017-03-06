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

import bioinf.mutations.Mutations
import Mutations._
import bioinf.Answer._
import bioinf.Output._
import scala.io.Source
val it = Source.fromFile("""C:\tmp\dataset_296_6.txt""").getLines()
val text1 = it.next()
val text2 = it.next()
val s = suffixTreeConstruction(text1 + "#" + text2 + "$")
val dc = deepestColoredEdges(s)
val lss = dc.map{printSuffix(s,_)}
val answer = lss.head
writeStringToFile(answer,"""c:\tmp\results.txt""")
