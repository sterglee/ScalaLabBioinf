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

import bioinf.hiddenmessages.HiddenMessages
import HiddenMessages._
import bioinf.Answer._
import scala.io.Source

val it = Source.fromFile("C:\\tmp\\dataset_4_5.txt").getLines()
val genome = it.next()
val args = it.next().split(" ").map(_.toInt)
assert(args.length == 3)
val k = args(0)
val L = args(1)
val t = args(2)

printSeq(clumpFindingNaive(genome,k,L,t))