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

import bioinf.UnitSpec
import HiddenMessages._

class HiddenMessagesSpec extends UnitSpec {

  "HiddenMessages" should "generate kmers" in {

    val kmers = allKmers(2)
    val all3mers = IndexedSeq[String](
      "AA","AC","AG","AT",
      "CA","CC","CG","CT",
      "GA","GC","GG","GT",
      "TA","TC","TG","TT"

    )

    kmers should be (all3mers)
  }
}
