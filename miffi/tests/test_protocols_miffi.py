# **************************************************************************
# *
# * Authors:     Daniel Marchan (da.marchan@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportMicrographs

from ..protocols import MiffiProtMicrographs
from ..protocols import protocol_miffi

import pickle
from pathlib import Path

class TestMiffi(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('xmipp_tutorial')
        cls.micsFn = cls.ds.getFile('micrographs/*.mrc')

        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.newProtocol(ProtImportMicrographs,
                                             filesPath=cls.micsFn,
                                             samplingRate=2.56)
        cls.launchProtocol(cls.protImportMics)

    def testMiffi(self):
        print(magentaStr("\n==> Testing miffi - categorize micrographs:"))
        protMiffi = self.newProtocol(MiffiProtMicrographs, inputSet=self.protImportMics.outputMicrographs)
        self.launchProtocol(protMiffi)
        micSet = getattr(protMiffi, 'outputMicrographsDiscarded', None)
        self.assertEqual(len(micSet), self._getNumberOfBad(protMiffi))

    def _getNumberOfBad(self, protMiffi):
        categorize_output_file = protocol_miffi.find_file_with_ending(protMiffi._getExtraPath('micBatch1'),'dict.pkl')
        all_bad = []

        with open(categorize_output_file, 'rb') as file:
            data = pickle.load(file)
            # Check for the presence of 'good_high_conf' and 'good_low_conf' and extract file names
            if 'bad_single' in data:
                all_bad.extend([Path(path).name for path in data['bad_single']])
            if 'bad_multiple' in data:
                all_bad.extend([Path(path).name for path in data['bad_multiple']])

        return len(all_bad)