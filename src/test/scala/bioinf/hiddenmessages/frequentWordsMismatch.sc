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
import bioinf.Input._
import bioinf.Output._
import scala.io.Source
val it = Source.fromFile("""C:\tmp\test.txt""").getLines()
val a = it.next()
val b = it.next().split(" ")
val answer = printSeq(frequentWordsMismatch(a,b.head.toInt,b.last.toInt))
writeStringToFile(answer,"""c:\tmp\results.txt""")