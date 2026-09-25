# ***************************************************************************
# *
# * Regression tests for MIFFI streaming/resume behavior.
# *
# ***************************************************************************

import unittest
from datetime import datetime
from unittest.mock import patch

from miffi.protocols.protocol_miffi import MiffiProtMicrographs
from miffi.protocols import protocol_miffi as miffi_module


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



class TestMiffiPendingResultsRegression(unittest.TestCase):
    def testResultsFinishedDuringOutputSnapshotAreNotLost(self):
        import copy as stdcopy
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict
        from unittest.mock import patch

        from miffi.protocols import protocol_miffi as miffi_module

        with tempfile.TemporaryDirectory() as tmp:
            firstPkl = os.path.join(tmp, "batch_1_dict.pkl")
            firstLog = os.path.join(tmp, "batch_1.log")
            latePkl = os.path.join(tmp, "batch_2_dict.pkl")
            lateLog = os.path.join(tmp, "batch_2.log")

            for path in (firstPkl, latePkl):
                with open(path, "wb") as handle:
                    pickle.dump({}, handle)

            for path in (firstLog, lateLog):
                with open(path, "w") as handle:
                    handle.write("")

            class _Summary:
                def __init__(self):
                    self.value = ""

                def set(self, value):
                    self.value = value

            class _Image:
                def __init__(self, objId):
                    self.objId = objId

                def clone(self):
                    return _Image(self.objId)

                def getFileName(self):
                    return os.path.join(tmp, f"mic_{self.objId}.mrc")

            class _Input:
                def getSize(self):
                    return 2

                def getItem(self, field, objId):
                    assert field == "id"
                    return _Image(objId)

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1]
                    self.outputCategorizeFiles = [firstPkl]
                    self.outputCategorizeLogFiles = [firstLog]
                    self.isStreamClosed = True
                    self.inputFn = "logical-input"
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = []
                    self.firstTime = {
                        miffi_module.OUTPUT: True,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.summaryVar = _Summary()
                    self.inputSet = object()

                def _getAllDoneIds(self):
                    return [], 0, [], []

                def _loadInputSet(self, inputFn):
                    return _Input()

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass

                def _getFirstJoinStep(self):
                    return None

                def _store(self):
                    pass

                def prepareBatch(self, newIds, counterBatch):
                    return os.path.join(tmp, f"batch_{counterBatch}")

                def runMiffiInference(self, batchDir, counterBatch):
                    return "inference.pkl"

                def runMiffiCategorize(self, counterBatch, inferenceFile):
                    return latePkl, lateLog

                def deleteBatch(self, batchDir):
                    pass

                def info(self, *args, **kwargs):
                    pass

                def delayRegister(self):
                    pass

            protocol = _Harness()
            snapshotStarted = threading.Event()
            workerDone = threading.Event()
            workerErrors = []
            realDeepcopy = stdcopy.deepcopy

            def _finishSecondBatch():
                try:
                    snapshotStarted.wait(timeout=2)
                    miffi_module.MiffiProtMicrographs.miffStep(
                        protocol,
                        [2],
                        2,
                    )
                except Exception as exc:
                    workerErrors.append(exc)
                finally:
                    workerDone.set()

            def _controlledDeepcopy(value, memo=None):
                if value is protocol.outputCategorizeFiles:
                    snapshot = realDeepcopy(value, memo)
                    snapshotStarted.set()

                    # Without synchronization miffStep can append right here,
                    # between the snapshot and the queue reset. With the
                    # protocol result lock held, the worker must wait until
                    # the snapshot/reset transaction is complete.
                    if not protocol._resultsLock.locked():
                        if not workerDone.wait(timeout=2):
                            raise AssertionError(
                                "The simulated MIFFI worker did not finish "
                                "inside the unsynchronized snapshot window."
                            )
                    return snapshot

                return realDeepcopy(value, memo)

            worker = threading.Thread(target=_finishSecondBatch)
            worker.start()

            try:
                with patch.object(
                    miffi_module.copy,
                    "deepcopy",
                    side_effect=_controlledDeepcopy,
                ):
                    miffi_module.MiffiProtMicrographs._checkNewOutput(
                        protocol
                    )
            finally:
                snapshotStarted.set()
                worker.join(timeout=2)

            self.assertFalse(
                worker.is_alive(),
                "The simulated MIFFI worker must not remain blocked.",
            )
            self.assertEqual([], workerErrors)

            self.assertIn(
                latePkl,
                protocol.outputCategorizeFiles,
                "A batch finishing while output results are snapshotted "
                "must remain pending for the next registration cycle.",
            )
            self.assertIn(
                lateLog,
                protocol.outputCategorizeLogFiles,
            )
            self.assertIn(2, protocol.processedIds)


class _ExistingOutputSet:
    def __init__(self):
        self.enableAppendCalls = 0
        self.copiedFrom = None

    def enableAppend(self):
        self.enableAppendCalls += 1

    def copyInfo(self, inputs):
        self.copiedFrom = inputs


class TestMiffiLoadOutputSetRegression(unittest.TestCase):
    def testLoadOutputSetReusesLogicalOutputWithoutBackingFile(self):
        # Regression test: an output that Scipion already knows about
        # (protocol.outputMicrographs) must be reused even when its backing
        # file was never materialized on disk yet. Falling through to "no
        # backing file -> build a fresh, empty Set" would silently discard
        # whatever was already appended to the real logical output.
        existingOutputSet = _ExistingOutputSet()
        inputs = object()

        class _Harness:
            def __init__(self):
                self.outputMicrographs = existingOutputSet
                self.inputSet = _Pointer(inputs)

            def _getPath(self, name):
                return "/tmp/" + name

        protocol = _Harness()

        with patch.object(miffi_module.os.path, "exists", return_value=False):
            outputSet = MiffiProtMicrographs._loadOutputSet(
                protocol, object, "micrographs.dat",
                outputName=miffi_module.OUTPUT,
            )

        self.assertIs(existingOutputSet, outputSet)
        self.assertEqual(1, existingOutputSet.enableAppendCalls)
        self.assertIs(inputs, existingOutputSet.copiedFrom)


if __name__ == "__main__":
    unittest.main()
