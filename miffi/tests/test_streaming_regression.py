# ***************************************************************************
# *
# * Regression tests for MIFFI streaming/resume behavior.
# *
# ***************************************************************************

import unittest
from datetime import datetime

from miffi.protocols.protocol_miffi import MiffiProtMicrographs


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _Pointer:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _LogicalInputSet:
    """Backend-agnostic logical Set used by the protocol pointer."""

    def __init__(self):
        self.closeCalls = 0
        self.loadCalls = 0
        self.loadPropertiesCalls = 0

    def getFileName(self):
        return "/tmp/input-micrographs.sqlite"

    def hasChangedSince(self, _lastCheck):
        return True

    def close(self):
        self.closeCalls += 1

    def load(self):
        self.loadCalls += 1
        return self

    def loadAllProperties(self):
        self.loadPropertiesCalls += 1

    def getIdSet(self):
        return {1, 2}

    def getSize(self):
        return 2

    def isStreamClosed(self):
        return False


class _StaleFileBackedSet:
    """Represents a storage snapshot that has not caught up yet."""

    def __init__(self, filename=None):
        self.filename = filename

    def loadAllProperties(self):
        pass

    def getIdSet(self):
        return {1}

    def getSize(self):
        return 1

    def isStreamClosed(self):
        return False

    def close(self):
        pass


class _StreamingHarness:
    def __init__(self):
        self.logicalInput = _LogicalInputSet()
        self.inputSet = _Pointer(self.logicalInput)
        self.inputFn = self.logicalInput.getFileName()
        self._inputClass = _StaleFileBackedSet
        self.streamingBatchSize = _Value(1)
        self.insertedIds = [1]
        self.isStreamClosed = False
        self.lastCheck = datetime.now()
        self.insertedBatches = []
        self.updateStepsCalls = 0

    def debug(self, _message):
        pass

    def _loadInputSet(self, inputFn):
        return MiffiProtMicrographs._loadInputSet(self, inputFn)

    def _getFirstJoinStep(self):
        return None

    def isContinued(self):
        return False

    def _insertNewImageSteps(self, newIds, batchSize):
        ids = list(newIds)
        self.insertedBatches.append(ids)
        self.insertedIds.extend(ids)
        return [101]

    def updateSteps(self):
        self.updateStepsCalls += 1


class TestMiffiStreamingRegression(unittest.TestCase):
    def testStreamingReloadsLogicalSetWhenStorageSnapshotDidNotChange(self):
        protocol = _StreamingHarness()

        MiffiProtMicrographs._checkNewInput(protocol)

        self.assertEqual(protocol.insertedBatches, [[2]])
        self.assertEqual(protocol.insertedIds, [1, 2])
        self.assertEqual(protocol.updateStepsCalls, 1)
        self.assertEqual(protocol.logicalInput.closeCalls, 2)
        self.assertEqual(protocol.logicalInput.loadCalls, 1)
        self.assertEqual(protocol.logicalInput.loadPropertiesCalls, 1)


class _FailingBatchHarness:
    def __init__(self):
        self.processedIds = []
        self.outputCategorizeFiles = []
        self.outputCategorizeLogFiles = []
        self.isStreamClosed = True
        self.deletedBatches = []
        self.messages = []

    def prepareBatch(self, newIds, counterBatch):
        return "/tmp/miffi-failing-batch"

    def runMiffiInference(self, batchDir, counterBatch):
        raise RuntimeError("simulated MIFFI failure")

    def runMiffiCategorize(self, counterBatch, inferenceFile):
        raise AssertionError("categorize must not run after inference failure")

    def info(self, message):
        self.messages.append(str(message))

    def deleteBatch(self, batchDir):
        self.deletedBatches.append(batchDir)

    def delayRegister(self):
        pass


class TestMiffiBatchFailureRegression(unittest.TestCase):
    def testFailedBatchIsNotReportedAsProcessed(self):
        protocol = _FailingBatchHarness()

        with self.assertRaises(RuntimeError):
            MiffiProtMicrographs.miffStep(protocol, [7, 8], 1)

        self.assertEqual(protocol.processedIds, [])
        self.assertEqual(protocol.outputCategorizeFiles, [])
        self.assertEqual(protocol.outputCategorizeLogFiles, [])


if __name__ == "__main__":
    unittest.main()
