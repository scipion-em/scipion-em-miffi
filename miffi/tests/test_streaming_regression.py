# ***************************************************************************
# *
# * Regression tests for MIFFI streaming/resume behavior.
# *
# ***************************************************************************

import os
import unittest
from datetime import datetime
from unittest.mock import patch

from pyworkflow.protocol import ProtStreamingBase
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
        self.newDeps = []

    def debug(self, _message):
        pass

    def _loadInputSet(self, inputFn):
        return MiffiProtMicrographs._loadInputSet(self, inputFn)

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

    def testStreamingPollDiscoversOnlyIdsAboveWatermark(self):
        class _IncrementalInputSet:
            def __init__(self):
                self.whereCalls = []
                self.closed = False

            def hasChangedSince(self, _lastCheck):
                return True

            def getIdSet(self):
                raise AssertionError(
                    "Normal MIFFI streaming polls must not full-scan all input ids."
                )

            def getUniqueValues(self, attributes, where=None):
                assert attributes == "id"
                self.whereCalls.append(where)
                if where == "id > 2":
                    return [3, 4]
                raise AssertionError(
                    "Expected an incremental id query above the watermark, got %r"
                    % where
                )

            def isStreamClosed(self):
                return False

            def close(self):
                self.closed = True

        class _BatchSize:
            def get(self):
                return 1

        class _Harness:
            def __init__(self):
                self._inputSet = _IncrementalInputSet()
                self.inputSet = _Pointer(self._inputSet)
                self.inputFn = "logical-input"
                self.insertedIds = [1, 2]
                self._inputWatermark = 2
                self._pendingInputIds = set()
                self.isStreamClosed = False
                self.streamingBatchSize = _BatchSize()
                self.newDeps = []
                self.scheduled = []
                self.lastCheck = None

            def debug(self, _message):
                pass

            def info(self, _message):
                pass

            def _loadInputSet(self, _inputFn):
                return self._inputSet

            def isContinued(self):
                return False

            def _insertNewImageSteps(self, newIds, batchSize):
                ids = list(newIds)
                self.scheduled.append((ids, batchSize))
                self.insertedIds.extend(ids)
                return [101]

            def updateSteps(self):
                pass

        protocol = _Harness()

        MiffiProtMicrographs._checkNewInput(protocol)

        self.assertEqual(
            ["id > 2"],
            protocol._inputSet.whereCalls,
            "MIFFI must query only ids above the current streaming watermark.",
        )
        self.assertEqual(
            [([3, 4], 1)],
            protocol.scheduled,
        )
        self.assertEqual(
            4,
            protocol._inputWatermark,
        )

    def testClosedStreamReconcilesIdBelowWatermark(self):
        class _ClosedInputSet:
            def __init__(self):
                self.whereCalls = []
                self.closed = False

            def getUniqueValues(self, attributes, where=None):
                assert attributes == "id"
                self.whereCalls.append(where)
                if where == "id > 10":
                    return []
                if where is None:
                    return [9, 10]
                raise AssertionError("Unexpected where clause: %r" % where)

            def getSize(self):
                return 2

            def isStreamClosed(self):
                return True

            def close(self):
                self.closed = True

        class _BatchSize:
            def get(self):
                return 1

        class _Harness:
            def __init__(self):
                self._inputSet = _ClosedInputSet()
                self.inputSet = _Pointer(self._inputSet)
                self.inputFn = "logical-input"
                self.insertedIds = [10]
                self._inputWatermark = 10
                self._pendingInputIds = set()
                self.isStreamClosed = False
                self.streamingBatchSize = _BatchSize()
                self.newDeps = []
                self.scheduled = []
                self.lastCheck = None

            def debug(self, _message):
                pass

            def info(self, _message):
                pass

            def _loadInputSet(self, _inputFn):
                return self._inputSet

            def isContinued(self):
                return False

            def _insertNewImageSteps(self, newIds, batchSize):
                ids = list(newIds)
                self.scheduled.append((ids, batchSize))
                self.insertedIds.extend(ids)
                return [202]

            def updateSteps(self):
                pass

        protocol = _Harness()

        MiffiProtMicrographs._checkNewInput(protocol)

        self.assertEqual(
            ["id > 10", None],
            protocol._inputSet.whereCalls,
            "A closed stream with fewer known ids than its logical size must "
            "perform one terminal full-id reconciliation.",
        )
        self.assertEqual(
            [([9], 1)],
            protocol.scheduled,
            "Terminal reconciliation must recover an id that became visible "
            "below the already advanced watermark.",
        )
        self.assertEqual(10, protocol._inputWatermark)

    def testReleasedLateVisibleIdBypassesNoChangeShortcut(self):
        protocol = _StreamingHarness()

        # Represents the state after mic 1 was already discovered but a
        # worker could not select it yet and released it from insertedIds.
        # Mic 2 remains scheduled, so insertedIds is still non-empty.
        protocol.insertedIds = [2]
        protocol.logicalInput.hasChangedSince = lambda _lastCheck: False

        MiffiProtMicrographs._checkNewInput(protocol)

        self.assertEqual(
            [[1]],
            protocol.insertedBatches,
            "A previously discovered mic released for a visibility retry "
            "must not be hidden by the no-change shortcut.",
        )


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

    def _prepareBatchWithIds(self, newIds, counterBatch):
        return self.prepareBatch(newIds, counterBatch), list(newIds)

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

                def getIdSet(self):
                    return {1, 2}

                def __contains__(self, objId):
                    return True

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


                def _store(self):
                    pass

                def prepareBatch(self, newIds, counterBatch):
                    return os.path.join(tmp, f"batch_{counterBatch}")

                def _prepareBatchWithIds(self, newIds, counterBatch):
                    return self.prepareBatch(newIds, counterBatch), list(newIds)

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

    def testOutputPreparationFailureKeepsResultFilesPending(self):
        import pickle
        import tempfile
        import threading

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_1_dict.pkl")
            logFile = os.path.join(tmp, "batch_1.log")

            with open(pklFile, "wb") as handle:
                pickle.dump({miffi_module.GOOD: ["mic_1.mrc"]}, handle)
            with open(logFile, "w") as handle:
                handle.write("")

            class _Input:
                def getSize(self):
                    return 1

                def getIdSet(self):
                    return {1}

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = True
                    self.inputFn = "logical-input"
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = []

                def _getAllDoneIds(self):
                    return [], 0, [], []

                def _loadInputSet(self, _inputFn):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    raise RuntimeError("simulated output preparation failure")

            protocol = _Harness()

            with self.assertRaises(RuntimeError):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                [pklFile],
                protocol.outputCategorizeFiles,
                "MIFFI classification results must remain pending when output preparation fails.",
            )
            self.assertEqual(
                [logFile],
                protocol.outputCategorizeLogFiles,
            )

    def testPartialOutputFailureKeepsAlreadyPersistedIdsDurable(self):
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_1_dict.pkl")
            logFile = os.path.join(tmp, "batch_1.log")

            with open(pklFile, "wb") as handle:
                pickle.dump(
                    {
                        miffi_module.GOOD: ["mic_1.mrc"],
                        miffi_module.BAD_FILM: ["mic_2.mrc"],
                    },
                    handle,
                )
            with open(logFile, "w") as handle:
                handle.write("")

            class _Image:
                def __init__(self, objId):
                    self.objId = objId

                def clone(self):
                    return _Image(self.objId)

                def getFileName(self):
                    return os.path.join(tmp, "mic_%d.mrc" % self.objId)

            class _Input:
                def __contains__(self, objId):
                    return objId in {1, 2}

                def getItem(self, field, objId):
                    assert field == "id"
                    return _Image(objId)

            class _Output:
                def __init__(self):
                    self.ids = set()

                def append(self, image):
                    self.ids.add(image.objId)

                def __contains__(self, objId):
                    return objId in self.ids

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1, 2]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = False
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = [miffi_module.BAD_FILM]
                    self.firstTime = {
                        miffi_module.OUTPUT: True,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.inputSet = object()
                    self.acceptedOutput = _Output()
                    self.discardedOutput = _Output()
                    # Simulate a normal long-running stream whose durable-id
                    # cache was already populated by an earlier poll.
                    self._persistedDoneIdsCache = (set(), set(), set())

                def _getAllDoneIds(self):
                    return MiffiProtMicrographs._getAllDoneIds(self)

                def _loadInputSet(self, _unused=None):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    if outputName == miffi_module.OUTPUT:
                        return self.acceptedOutput
                    return self.discardedOutput

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    if outputName == miffi_module.OUTPUT:
                        setattr(self, outputName, outputSet)
                        return
                    raise RuntimeError(
                        "simulated discarded-output persistence failure"
                    )

                def _defineSourceRelation(self, *_args):
                    pass

                def error(self, _message):
                    pass

            protocol = _Harness()

            with self.assertRaises(RuntimeError):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                [2],
                protocol.processedIds,
                "An accepted id that was already durably published before "
                "a later discarded-output failure must leave the pending "
                "processed queue.",
            )
            self.assertEqual(
                {1},
                protocol._persistedDoneIdsCache[0],
                "Partial publication must update the durable-id cache with "
                "the output that succeeded before the later failure.",
            )
            self.assertEqual([pklFile], protocol.outputCategorizeFiles)
            self.assertEqual([logFile], protocol.outputCategorizeLogFiles)


    def testPartialRetryUpdatesOnlyStillPendingOutput(self):
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_1_dict.pkl")
            logFile = os.path.join(tmp, "batch_1.log")

            with open(pklFile, "wb") as handle:
                pickle.dump(
                    {
                        miffi_module.GOOD: ["mic_1.mrc"],
                        miffi_module.BAD_FILM: ["mic_2.mrc"],
                    },
                    handle,
                )
            with open(logFile, "w") as handle:
                handle.write("")

            class _Summary:
                def set(self, _value):
                    pass

            class _Image:
                def __init__(self, objId):
                    self.objId = objId

                def clone(self):
                    return _Image(self.objId)

                def getFileName(self):
                    return os.path.join(tmp, "mic_%d.mrc" % self.objId)

            class _Input:
                def __contains__(self, objId):
                    return objId in {1, 2}

                def getItem(self, field, objId):
                    assert field == "id"
                    return _Image(objId)

            class _Output:
                def __init__(self):
                    self.ids = set()

                def append(self, image):
                    self.ids.add(image.objId)

                def __contains__(self, objId):
                    return objId in self.ids

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1, 2]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = False
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = [miffi_module.BAD_FILM]
                    self.firstTime = {
                        miffi_module.OUTPUT: True,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.summaryVar = _Summary()
                    self.inputSet = object()
                    self.acceptedOutput = _Output()
                    self.discardedOutput = _Output()
                    self._persistedDoneIdsCache = (set(), set(), set())
                    self.acceptedUpdateCalls = 0
                    self.discardedUpdateCalls = 0
                    self.failDiscardedOnce = True

                def _getAllDoneIds(self):
                    return MiffiProtMicrographs._getAllDoneIds(self)

                def _loadInputSet(self, _unused=None):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    if outputName == miffi_module.OUTPUT:
                        return self.acceptedOutput
                    return self.discardedOutput

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    if outputName == miffi_module.OUTPUT:
                        self.acceptedUpdateCalls += 1
                        setattr(self, outputName, outputSet)
                        return

                    self.discardedUpdateCalls += 1
                    if self.failDiscardedOnce:
                        self.failDiscardedOnce = False
                        raise RuntimeError(
                            "simulated discarded-output persistence failure"
                        )
                    setattr(self, outputName, outputSet)

                def _defineSourceRelation(self, *_args):
                    pass

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass

                def _store(self):
                    pass

                def error(self, _message):
                    pass

            protocol = _Harness()

            with self.assertRaises(RuntimeError):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual([2], protocol.processedIds)
            self.assertEqual({1}, protocol._persistedDoneIdsCache[0])

            MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                1,
                protocol.acceptedUpdateCalls,
                "Retrying a partially persisted MIFFI batch must not update "
                "an output whose candidate ids are already durable.",
            )
            self.assertEqual(2, protocol.discardedUpdateCalls)
            self.assertEqual([], protocol.processedIds)
            self.assertEqual({1, 2}, protocol._persistedDoneIdsCache[0])


    def testOutputUpdateFailureKeepsResultFilesPending(self):
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_1_dict.pkl")
            logFile = os.path.join(tmp, "batch_1.log")

            with open(pklFile, "wb") as handle:
                pickle.dump({miffi_module.GOOD: ["mic_1.mrc"]}, handle)
            with open(logFile, "w") as handle:
                handle.write("")

            class _Image:
                def clone(self):
                    return self

                def getFileName(self):
                    return os.path.join(tmp, "mic_1.mrc")

            class _Input:
                def getSize(self):
                    return 1

                def getIdSet(self):
                    return {1}

                def __contains__(self, objId):
                    return objId == 1

                def getItem(self, field, objId):
                    assert field == "id"
                    assert objId == 1
                    return _Image()

            class _Output:
                def append(self, _image):
                    pass

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
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
                    self.inputSet = object()

                def _getAllDoneIds(self):
                    return [], 0, [], []

                def _loadInputSet(self, _inputFn):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    return _Output()

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    raise RuntimeError("simulated output persistence failure")

                def error(self, _message):
                    pass

            protocol = _Harness()

            with self.assertRaises(RuntimeError):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                [pklFile],
                protocol.outputCategorizeFiles,
                "MIFFI classification results must remain pending when output persistence fails.",
            )
            self.assertEqual(
                [logFile],
                protocol.outputCategorizeLogFiles,
            )

    def testCheckNewOutputSkipsImageNotYetVisibleWithoutCrashing(self):
        # Regression test: Set.getItem raises (UnboundLocalError) rather
        # than returning None for a row it cannot find. An imageId
        # already in processedIds (because its batch step already
        # succeeded) may still momentarily fail to be selectable from a
        # freshly-reloaded input Set - it must be skipped (and logged)
        # instead of crashing the whole protocol, while the rest of the
        # batch is still processed normally.
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            firstPkl = os.path.join(tmp, "batch_1_dict.pkl")
            firstLog = os.path.join(tmp, "batch_1.log")
            with open(firstPkl, "wb") as handle:
                pickle.dump({miffi_module.GOOD: ["mic_1.mrc"]}, handle)
            with open(firstLog, "w") as handle:
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
                def __init__(self, visibleIds):
                    self._visibleIds = visibleIds

                def getSize(self):
                    return len(self._visibleIds)

                def getIdSet(self):
                    return set(self._visibleIds)

                def __contains__(self, objId):
                    return objId in self._visibleIds

                def getItem(self, field, objId):
                    assert field == "id"
                    return _Image(objId)

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1, 2]  # id 2 will not be visible
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
                    self.errors = []
                    self.appended = []
                    self._inputClass = _Image
                    self._baseName = "micrographs.sqlite"

                def _getAllDoneIds(self):
                    return [], 0, [], []

                def _loadInputSet(self, inputFn):
                    return _Input(visibleIds={1})

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass


                def _store(self):
                    pass

                def error(self, msg):
                    self.errors.append(msg)

                def info(self, *args, **kwargs):
                    pass

                def _loadOutputSet(self, outputName, suffix=""):
                    return self

                def append(self, image):
                    self.appended.append(image.objId)

                def __contains__(self, objId):
                    return objId in self.appended

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    pass

                def _defineSourceRelation(self, *args, **kwargs):
                    pass

            protocol = _Harness()

            miffi_module.MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                1,
                len(protocol.errors),
                "The invisible image must be logged, not silently ignored.",
            )
            self.assertEqual(
                [1],
                protocol.appended,
                "The still-visible image must still be processed normally.",
            )


class _ExistingOutputSet:
    def __init__(self):
        self.enableAppendCalls = 0
        self.loadAllPropertiesCalls = 0
        self.copiedFrom = None

    def loadAllProperties(self):
        self.loadAllPropertiesCalls += 1

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

        outputSet = MiffiProtMicrographs._loadOutputSet(protocol, miffi_module.OUTPUT)

        self.assertIs(existingOutputSet, outputSet)
        self.assertEqual(1, existingOutputSet.loadAllPropertiesCalls)
        self.assertEqual(1, existingOutputSet.enableAppendCalls)
        self.assertIs(inputs, existingOutputSet.copiedFrom)


class _FakeMic:
    def __init__(self, objId, tmp):
        self._objId = objId
        self._tmp = tmp

    def clone(self):
        return _FakeMic(self._objId, self._tmp)

    def getFileName(self):
        micFn = os.path.join(self._tmp, "mic_%d.mrc" % self._objId)
        if not os.path.exists(micFn):
            with open(micFn, "w"):
                pass
        return micFn


class _FakeMicSet:
    def __init__(self, ids, tmp):
        self._ids = set(ids)
        self._tmp = tmp

    def __contains__(self, micId):
        return micId in self._ids

    def getItem(self, field, micId):
        assert field == "id"
        # Real Set.getItem raises (UnboundLocalError) rather than
        # returning None for a row it cannot find - match that here so a
        # missing membership guard is caught by these tests.
        if micId not in self._ids:
            raise UnboundLocalError("row not found for id %r" % micId)
        return _FakeMic(micId, self._tmp)


class _PrepareBatchHarness:
    MIC_VISIBILITY_MAX_ATTEMPTS = 3
    MIC_VISIBILITY_RETRY_DELAY = 0

    def __init__(self, tmp, inputSet):
        self._tmp = tmp
        self._inputSet = inputSet
        self.inputFn = "/tmp/fake_mics.sqlite"
        self.errors = []

    def _getTmpPath(self, name):
        return os.path.join(self._tmp, name)

    def _loadInputSet(self, inputFn):
        return self._inputSet

    def _prepareBatchWithIds(self, newIds, counterBatch):
        return MiffiProtMicrographs._prepareBatchWithIds(
            self, newIds, counterBatch
        )

    def error(self, msg):
        self.errors.append(msg)


class TestMiffiPrepareBatchVisibilityRegression(unittest.TestCase):
    """Regression tests for Set.getItem visibility handling in prepareBatch."""

    def testPrepareBatchRetriesUntilMicBecomesVisible(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            inputSet = _FakeMicSet(ids=[], tmp=tmp)  # mic 1 starts invisible
            harness = _PrepareBatchHarness(tmp, inputSet)

            def becomeVisibleOnSleep(_delay):
                inputSet._ids = {1}

            with patch(
                "miffi.protocols.protocol_miffi.time.sleep",
                side_effect=becomeVisibleOnSleep,
            ):
                batchDir = MiffiProtMicrographs.prepareBatch(harness, [1], 1)

            self.assertEqual([], harness.errors)
            copiedFiles = os.listdir(batchDir)
            self.assertEqual(["1__mic_1.mrc"], copiedFiles)

    def testPrepareBatchExcludesMicThatNeverBecomesVisible(self):
        # Regression test: Set.getItem raises (UnboundLocalError) rather
        # than returning None for a row it cannot find. A micId that
        # never becomes visible must be excluded from the batch (and
        # logged) instead of crashing the whole one-shot worker step.
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            inputSet = _FakeMicSet(ids=[], tmp=tmp)  # never becomes visible
            harness = _PrepareBatchHarness(tmp, inputSet)

            with patch(
                "miffi.protocols.protocol_miffi.time.sleep",
                return_value=None,
            ):
                batchDir = MiffiProtMicrographs.prepareBatch(harness, [1], 1)

            self.assertEqual(1, len(harness.errors))
            self.assertEqual([], os.listdir(batchDir))


if __name__ == "__main__":
    unittest.main()


class _RefreshRequiredMiffiOutput:
    def __init__(self, ids):
        self.ids = set(ids)
        self.loaded = False
        self.appendEnabled = False
        self.copiedFrom = None

    def loadAllProperties(self):
        self.loaded = True

    def getSize(self):
        if not self.loaded:
            raise AssertionError("Persisted MIFFI output must be refreshed before reading its size.")
        return len(self.ids)

    def getIdSet(self):
        if not self.loaded:
            raise AssertionError("Persisted MIFFI output must be refreshed before reading its ids.")
        return set(self.ids)

    def enableAppend(self):
        if not self.loaded:
            raise AssertionError("Persisted MIFFI output must be refreshed before enableAppend().")
        self.appendEnabled = True

    def copyInfo(self, inputs):
        self.copiedFrom = inputs


class TestMiffiLogicalInputIdentityRegression(unittest.TestCase):
    def testInitializeParamsDoesNotReadInputBackingFileName(self):
        class _Input:
            def isStreamClosed(self):
                return False

            def getFileName(self):
                raise AssertionError(
                    "MIFFI streaming initialization must not depend on the "
                    "input Set backing filename."
                )

            def getClass(self):
                raise AssertionError(
                    "Unused input implementation class must not be captured "
                    "during streaming initialization."
                )

            def getClassName(self):
                raise AssertionError(
                    "Unused input implementation type must not be captured "
                    "during streaming initialization."
                )

        class _Harness:
            def __init__(self):
                self.inputSet = _Pointer(_Input())

            def _getDefinedLabels(self):
                return [], []

            def _outputNeedsSourceRelation(self, outputName):
                return MiffiProtMicrographs._outputNeedsSourceRelation(
                    self, outputName
                )

        protocol = _Harness()

        MiffiProtMicrographs.initializeParams(protocol)

        self.assertFalse(protocol.isStreamClosed)
        self.assertFalse(
            hasattr(protocol, "inputFn"),
            "MIFFI must not persist a backing filename as logical input identity.",
        )
        self.assertFalse(hasattr(protocol, "_inputClass"))
        self.assertFalse(hasattr(protocol, "_inputType"))


class TestMiffiSourceRelationFailureRegression(unittest.TestCase):
    def testPersistedOutputRetriesMissingSourceRelationWithoutRepublishing(self):
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_dict.pkl")
            logFile = os.path.join(tmp, "batch.log")

            with open(pklFile, "wb") as handle:
                pickle.dump(
                    {
                        miffi_module.GOOD: ["mic_1.mrc"],
                    },
                    handle,
                )
            with open(logFile, "w") as handle:
                handle.write("")

            class _Image:
                def clone(self):
                    return self

                def getFileName(self):
                    return os.path.join(tmp, "mic_1.mrc")

            class _Input:
                def __contains__(self, objId):
                    return objId == 1

                def getItem(self, field, objId):
                    assert field == "id"
                    assert objId == 1
                    return _Image()

            class _Output:
                def __init__(self):
                    self.ids = set()

                def append(self, image):
                    self.ids.add(1)

                def __contains__(self, objId):
                    return objId in self.ids

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = False
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = []
                    self.firstTime = {
                        miffi_module.OUTPUT: True,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.inputSet = object()
                    self.finished = False
                    self.output = _Output()
                    self._persistedDoneIdsCache = (set(), set(), set())
                    self.outputUpdateCalls = 0
                    self.relationCalls = 0
                    self.storeCalls = 0

                def _getAllDoneIds(self):
                    return MiffiProtMicrographs._getAllDoneIds(self)

                def _retryPendingSourceRelations(self):
                    return MiffiProtMicrographs._retryPendingSourceRelations(
                        self
                    )

                def _loadInputSet(self, _unused=None):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    assert outputName == miffi_module.OUTPUT
                    return self.output

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    self.outputUpdateCalls += 1
                    setattr(self, outputName, outputSet)

                def _defineSourceRelation(self, source, target):
                    self.relationCalls += 1
                    if self.relationCalls == 1:
                        raise RuntimeError(
                            "simulated source-relation persistence failure"
                        )

                def _store(self):
                    self.storeCalls += 1

                def error(self, _message):
                    pass

            protocol = _Harness()

            with self.assertRaises(RuntimeError):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                {1},
                protocol._persistedDoneIdsCache[0],
                "The successfully persisted output id must remain durable even "
                "when relation creation fails afterwards.",
            )
            self.assertEqual(
                [],
                protocol.processedIds,
                "A durably published id must not be republished just because "
                "source-relation creation failed.",
            )
            self.assertEqual(1, protocol.outputUpdateCalls)
            self.assertEqual(1, protocol.relationCalls)

            MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                1,
                protocol.outputUpdateCalls,
                "Retrying a missing source relation must not republish an "
                "already durable output candidate.",
            )
            self.assertEqual(
                2,
                protocol.relationCalls,
                "A source relation that failed before being inserted must be "
                "retried even when there are no new output candidates.",
            )


    def testPersistedSourceRelationIsNotDuplicatedAfterAmbiguousFailure(self):
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_dict.pkl")
            logFile = os.path.join(tmp, "batch.log")

            with open(pklFile, "wb") as handle:
                pickle.dump(
                    {
                        miffi_module.GOOD: ["1__mic_1.mrc"],
                    },
                    handle,
                )
            with open(logFile, "w") as handle:
                handle.write("")

            class _Image:
                def clone(self):
                    return self

                def getFileName(self):
                    return os.path.join(tmp, "mic_1.mrc")

            class _Input:
                def __contains__(self, objId):
                    return objId == 1

                def getItem(self, field, objId):
                    assert field == "id"
                    assert objId == 1
                    return _Image()

                def getObjId(self):
                    return 201

            class _Pointer:
                def __init__(self, obj):
                    self.obj = obj

                def get(self):
                    return self.obj

            class _Output:
                def __init__(self):
                    self.ids = set()

                def append(self, image):
                    self.ids.add(1)

                def __contains__(self, objId):
                    return objId in self.ids

                def getObjId(self):
                    return 301

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = False
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = []
                    self.firstTime = {
                        miffi_module.OUTPUT: True,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.logicalInput = _Input()
                    self.inputSet = _Pointer(self.logicalInput)
                    self.finished = False
                    self.output = _Output()
                    self._persistedDoneIdsCache = (set(), set(), set())
                    self.outputUpdateCalls = 0
                    self.relationCalls = 0
                    self.relationRows = []
                    self.storeCalls = 0

                def _getAllDoneIds(self):
                    return MiffiProtMicrographs._getAllDoneIds(self)

                def _outputNeedsSourceRelation(self, outputName):
                    return MiffiProtMicrographs._outputNeedsSourceRelation(
                        self, outputName
                    )

                def _retryPendingSourceRelations(self):
                    return MiffiProtMicrographs._retryPendingSourceRelations(
                        self
                    )

                def _loadInputSet(self, _unused=None):
                    return self.logicalInput

                def _loadOutputSet(self, outputName, suffix=""):
                    assert outputName == miffi_module.OUTPUT
                    return self.output

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    self.outputUpdateCalls += 1
                    setattr(self, outputName, outputSet)

                def getRelations(self):
                    return list(self.relationRows)

                def _defineSourceRelation(self, source, target):
                    self.relationCalls += 1
                    self.relationRows.append({
                        "name": miffi_module.RELATION_SOURCE,
                        "object_parent_id": 201,
                        "object_child_id": 301,
                    })
                    if self.relationCalls == 1:
                        raise RuntimeError(
                            "simulated failure after relation persistence"
                        )

                def _store(self):
                    self.storeCalls += 1

                def error(self, _message):
                    pass

            protocol = _Harness()

            with self.assertRaises(RuntimeError):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(1, protocol.relationCalls)
            self.assertEqual(1, len(protocol.relationRows))
            self.assertTrue(protocol.firstTime[miffi_module.OUTPUT])

            MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                1,
                protocol.relationCalls,
                "A datasource relation that is already durable must not be "
                "inserted again after an ambiguous persistence failure.",
            )
            self.assertEqual(
                1,
                len(protocol.relationRows),
                "Native Scipion relation storage has no uniqueness guard, so "
                "an ambiguous retry must not create a duplicate relation row.",
            )
            self.assertFalse(protocol.firstTime[miffi_module.OUTPUT])


class TestMiffiResumeOutputRelationRegression(unittest.TestCase):
    def testInitializeParamsDoesNotTreatExistingOutputAsFirstPublication(self):
        class _Input:
            def isStreamClosed(self):
                return False

        class _Harness:
            def __init__(self):
                self.inputSet = _Pointer(_Input())
                # Simulate Resume: the accepted logical output already exists
                # from the previous execution, while discarded does not.
                self.outputMicrographs = object()

            def _getDefinedLabels(self):
                return [], []

            def _outputNeedsSourceRelation(self, outputName):
                return MiffiProtMicrographs._outputNeedsSourceRelation(
                    self, outputName
                )

        protocol = _Harness()

        MiffiProtMicrographs.initializeParams(protocol)

        self.assertFalse(
            protocol.firstTime[miffi_module.OUTPUT],
            "An existing logical MIFFI output must not be treated as a first "
            "publication after Resume, otherwise the source relation is "
            "inserted again.",
        )
        self.assertTrue(
            protocol.firstTime[miffi_module.OUTPUT_DISCARDED],
            "An output that does not exist yet must still create its source "
            "relation on first publication.",
        )


class TestMiffiPersistedOutputRefresh(unittest.TestCase):
    def testDoneIdsRefreshPersistedAcceptedAndDiscardedOutputs(self):
        class Harness:
            def __init__(self):
                self.outputMicrographs = _RefreshRequiredMiffiOutput([1, 2])
                self.outputMicrographsDiscarded = _RefreshRequiredMiffiOutput([3])

        protocol = Harness()
        doneIds, sizeOutput, acceptedIds, discardedIds = MiffiProtMicrographs._getAllDoneIds(protocol)

        self.assertEqual({1, 2, 3}, set(doneIds))
        self.assertEqual(3, sizeOutput)
        self.assertEqual({1, 2}, set(acceptedIds))
        self.assertEqual({3}, set(discardedIds))
        self.assertTrue(protocol.outputMicrographs.loaded)
        self.assertTrue(protocol.outputMicrographsDiscarded.loaded)

    def testIdleOutputPollsRestorePersistedIdsOnlyOnce(self):
        import threading

        class _Output:
            def __init__(self):
                self.getIdSetCalls = 0

            def loadAllProperties(self):
                pass

            def getSize(self):
                return 1

            def getIdSet(self):
                self.getIdSetCalls += 1
                return {1}

        class _Input:
            def getIdSet(self):
                return {1, 2}

        class _Harness:
            def __init__(self):
                self._resultsLock = threading.Lock()
                self.outputMicrographs = _Output()
                self.processedIds = []
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self.isStreamClosed = False
                self.inputFn = "logical-input"
                self.finished = False
                self.firstTime = {
                    miffi_module.OUTPUT: False,
                    miffi_module.OUTPUT_DISCARDED: False,
                }
                self.storeCalls = 0

            def _getAllDoneIds(self):
                return MiffiProtMicrographs._getAllDoneIds(self)

            def _loadInputSet(self, _inputFn):
                return _Input()

            def _retryPendingSourceRelations(self):
                return MiffiProtMicrographs._retryPendingSourceRelations(
                    self
                )

            def _store(self):
                self.storeCalls += 1

        protocol = _Harness()

        MiffiProtMicrographs._checkNewOutput(protocol)
        MiffiProtMicrographs._checkNewOutput(protocol)

        self.assertEqual(
            1,
            protocol.outputMicrographs.getIdSetCalls,
            "Idle MIFFI output polls must restore all persisted output ids "
            "only once instead of full-scanning the output Set every poll.",
        )
        self.assertEqual(2, protocol.storeCalls)

    def testOpenIdleOutputPollDoesNotFullScanInputIds(self):
        import threading

        class _Input:
            def getIdSet(self):
                raise AssertionError(
                    "Open MIFFI output polls must not full-scan all input ids."
                )

        class _Harness:
            def __init__(self):
                self._resultsLock = threading.Lock()
                self.processedIds = []
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self.isStreamClosed = False
                self.inputFn = "logical-input"
                self.finished = False
                self.firstTime = {
                    miffi_module.OUTPUT: False,
                    miffi_module.OUTPUT_DISCARDED: False,
                }
                self.storeCalls = 0

            def _getAllDoneIds(self):
                return [], 0, [], []

            def _loadInputSet(self, _inputFn):
                return _Input()

            def _retryPendingSourceRelations(self):
                return MiffiProtMicrographs._retryPendingSourceRelations(
                    self
                )

            def _store(self):
                self.storeCalls += 1

        protocol = _Harness()

        MiffiProtMicrographs._checkNewOutput(protocol)

        self.assertFalse(
            protocol.finished,
            "An open stream cannot become terminal during an idle output poll.",
        )
        self.assertEqual(1, protocol.storeCalls)

    def testClosedStreamReconcilesPersistedOutputsBeforeCompletion(self):
        import threading

        class _Output:
            def __init__(self):
                self.getIdSetCalls = 0

            def loadAllProperties(self):
                pass

            def getIdSet(self):
                self.getIdSetCalls += 1
                return {1}

        class _Input:
            def getIdSet(self):
                return {1}

        class _Harness:
            def __init__(self):
                self._resultsLock = threading.Lock()
                self.outputMicrographs = _Output()
                self._persistedDoneIdsCache = (set(), set(), set())
                self.processedIds = []
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self.isStreamClosed = True
                self.inputFn = "logical-input"
                self.finished = False
                self.firstTime = {
                    miffi_module.OUTPUT: False,
                    miffi_module.OUTPUT_DISCARDED: False,
                }
                self.storeCalls = 0

            def _getAllDoneIds(self):
                return MiffiProtMicrographs._getAllDoneIds(self)

            def _loadInputSet(self, _inputFn):
                return _Input()

            def _retryPendingSourceRelations(self):
                return MiffiProtMicrographs._retryPendingSourceRelations(
                    self
                )

            def _store(self):
                self.storeCalls += 1

        protocol = _Harness()

        MiffiProtMicrographs._checkNewOutput(protocol)

        self.assertEqual(
            1,
            protocol.outputMicrographs.getIdSetCalls,
            "A closed MIFFI stream must reconcile persisted output ids "
            "instead of trusting a possibly stale in-memory cache.",
        )
        self.assertTrue(
            protocol.finished,
            "Terminal completion must use the reconciled durable output ids.",
        )
        self.assertEqual(
            {1},
            protocol._persistedDoneIdsCache[0],
        )
        self.assertEqual(1, protocol.storeCalls)

    def testDurablyPublishedIdsAreRemovedFromProcessedQueue(self):
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_dict.pkl")
            logFile = os.path.join(tmp, "batch.log")

            with open(pklFile, "wb") as handle:
                pickle.dump(
                    {miffi_module.GOOD: ["mic_2.mrc"]},
                    handle,
                )
            with open(logFile, "w") as handle:
                handle.write("")

            class _Summary:
                def set(self, _value):
                    pass

            class _Image:
                def __init__(self, objId):
                    self.objId = objId

                def clone(self):
                    return _Image(self.objId)

                def getFileName(self):
                    return os.path.join(tmp, "mic_%d.mrc" % self.objId)

            class _Input:
                def __contains__(self, objId):
                    return objId == 2

                def getItem(self, field, objId):
                    assert field == "id"
                    return _Image(objId)

            class _Output:
                STREAM_OPEN = miffi_module.Set.STREAM_OPEN

                def __init__(self):
                    self.ids = set()

                def append(self, image):
                    self.ids.add(image.objId)

                def __contains__(self, objId):
                    return objId in self.ids

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [2]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = False
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
                    self.finished = False
                    self.output = _Output()

                def _getAllDoneIds(self):
                    return [], 0, [], []

                def _loadInputSet(self, _inputFn):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    return self.output

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    setattr(self, outputName, outputSet)

                def _defineSourceRelation(self, *_args):
                    pass

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass

                def _store(self):
                    pass

                def error(self, _message):
                    pass

            protocol = _Harness()

            with patch.object(
                miffi_module,
                "setLabel",
                side_effect=lambda *_args, **_kwargs: None,
            ), patch.object(
                miffi_module,
                "populate_and_update_categories",
                return_value={},
            ):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                [],
                protocol.processedIds,
                "Durably published MIFFI ids must leave processedIds so "
                "future output polls do not rescan the full processing history.",
            )

    def testPublishedBatchVerifiesOnlyCandidateOutputIds(self):
        import os
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch_dict.pkl")
            logFile = os.path.join(tmp, "batch.log")

            with open(pklFile, "wb") as handle:
                pickle.dump(
                    {miffi_module.GOOD: ["mic_2.mrc"]},
                    handle,
                )
            with open(logFile, "w") as handle:
                handle.write("")

            class _Summary:
                def set(self, _value):
                    pass

            class _Image:
                def __init__(self, objId):
                    self.objId = objId

                def clone(self):
                    return _Image(self.objId)

                def getFileName(self):
                    return os.path.join(
                        tmp,
                        "mic_%d.mrc" % self.objId,
                    )

            class _Input:
                def getIdSet(self):
                    return {1, 2}

                def __contains__(self, objId):
                    return objId == 2

                def getItem(self, field, objId):
                    assert field == "id"
                    assert objId == 2
                    return _Image(objId)

            class _Output:
                STREAM_OPEN = miffi_module.Set.STREAM_OPEN

                def __init__(self):
                    self.ids = {1}
                    self.getIdSetCalls = 0
                    self.candidateChecks = []

                def loadAllProperties(self):
                    pass

                def getIdSet(self):
                    self.getIdSetCalls += 1
                    return set(self.ids)

                def enableAppend(self):
                    pass

                def copyInfo(self, _inputSet):
                    pass

                def append(self, image):
                    self.ids.add(image.objId)

                def __contains__(self, objId):
                    self.candidateChecks.append(objId)
                    return objId in self.ids

            class _Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.outputMicrographs = _Output()
                    self.processedIds = [2]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = False
                    self.inputFn = "logical-input"
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = []
                    self.firstTime = {
                        miffi_module.OUTPUT: False,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.summaryVar = _Summary()
                    self.inputSet = _Pointer(object())
                    self.finished = False

                def _getAllDoneIds(self):
                    return MiffiProtMicrographs._getAllDoneIds(self)

                def _loadInputSet(self, _inputFn):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    return MiffiProtMicrographs._loadOutputSet(
                        self,
                        outputName,
                        suffix,
                    )

                def _updateOutputSet(
                        self,
                        outputName,
                        outputSet,
                        _streamMode,
                ):
                    setattr(self, outputName, outputSet)

                def _defineSourceRelation(self, *_args):
                    pass

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass

                def _store(self):
                    pass

                def error(self, _message):
                    pass

            protocol = _Harness()

            with patch.object(
                miffi_module,
                "setLabel",
                side_effect=lambda *_args, **_kwargs: None,
            ), patch.object(
                miffi_module,
                "populate_and_update_categories",
                return_value={},
            ):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertEqual(
                1,
                protocol.outputMicrographs.getIdSetCalls,
                "Publishing a MIFFI batch must not invalidate the persisted "
                "output cache and full-scan the complete output Set again.",
            )
            self.assertEqual(
                [2],
                protocol.outputMicrographs.candidateChecks,
                "Durability after publishing must be verified with a "
                "point lookup for the newly published candidate id.",
            )

    def testLoadOutputSetRefreshesExistingLogicalOutputBeforeAppend(self):
        existingOutput = _RefreshRequiredMiffiOutput([1])
        inputs = object()

        class Harness:
            def __init__(self):
                self.outputMicrographs = existingOutput
                self.inputSet = _Pointer(inputs)

        protocol = Harness()
        outputSet = MiffiProtMicrographs._loadOutputSet(protocol, miffi_module.OUTPUT)

        self.assertIs(existingOutput, outputSet)
        self.assertTrue(existingOutput.loaded)
        self.assertTrue(existingOutput.appendEnabled)
        self.assertIs(inputs, existingOutput.copiedFrom)


class TestMiffiTerminalPersistenceRegression(unittest.TestCase):
    def testProcessedButUnpersistedIdDoesNotCloseStreaming(self):
        import threading
        from collections import defaultdict

        class _Summary:
            def set(self, value):
                pass

        class _Input:
            def getSize(self):
                return 2

            def getIdSet(self):
                return {1, 2}

            def __contains__(self, objId):
                return objId == 1

        class _Harness:
            def __init__(self):
                self._resultsLock = threading.Lock()
                self.processedIds = [2]
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self.isStreamClosed = True
                self.inputFn = "logical-input"
                self.acceptedLabels = [miffi_module.GOOD]
                self.rejectedLabels = []
                self.firstTime = {miffi_module.OUTPUT: True, miffi_module.OUTPUT_DISCARDED: True}
                self.labelHistory = defaultdict(list)
                self.timeHistory = []
                self.outputLog = {}
                self.summaryVar = _Summary()
                self.inputSet = object()
                self.finished = False
                self.errors = []

            def _getAllDoneIds(self):
                return [1], 1, [1], []

            def _loadInputSet(self, inputFn):
                return _Input()

            def _plotMiffiLabelHistogram(self):
                pass

            def _plotMiffiTimeEvolution(self):
                pass

            def _store(self):
                pass

            def error(self, message):
                self.errors.append(message)

        protocol = _Harness()

        with patch.object(miffi_module, "populate_and_update_categories", return_value={}):
            MiffiProtMicrographs._checkNewOutput(protocol)

        # finished is now the only completion signal: the streaming
        # generator's while-loop reads it directly, so it staying False is
        # exactly what keeps the protocol alive for another poll.
        self.assertFalse(protocol.finished, "A processed id must not make the protocol terminal until it is durably present in an output Set.")


class TestMiffiTerminalCompletionRegression(unittest.TestCase):
    def testClosedStreamRequiresExactInputIdsToBePersisted(self):
        import threading

        class _Input:
            def getSize(self):
                return 2

            def getIdSet(self):
                return {1, 2}

        class _Harness:
            def __init__(self):
                self._resultsLock = threading.Lock()
                self.processedIds = []
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self.isStreamClosed = True
                self.inputFn = "logical-input"
                self.finished = False
                self.firstTime = {
                    miffi_module.OUTPUT: False,
                    miffi_module.OUTPUT_DISCARDED: False,
                }

            def _getAllDoneIds(self):
                # Same cardinality as the input, but not the same logical IDs.
                return [1, 99], 2, [1, 99], []

            def _loadInputSet(self, _inputFn):
                return _Input()

            def _retryPendingSourceRelations(self):
                return MiffiProtMicrographs._retryPendingSourceRelations(
                    self
                )

            def _store(self):
                pass

        protocol = _Harness()

        MiffiProtMicrographs._checkNewOutput(protocol)

        self.assertFalse(
            protocol.finished,
            "A closed MIFFI stream must not become terminal merely because "
            "the number of persisted output IDs matches the input size; "
            "the persisted IDs must exactly match the logical input IDs.",
        )


class TestMiffiUnclassifiedResultRegression(unittest.TestCase):
    def testProcessedMicWithoutRecognizedCategoryIsReleasedForRetry(self):
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(
                tmp,
                "batch.pkl",
            )
            logFile = os.path.join(
                tmp,
                "batch.log",
            )

            with open(
                    pklFile,
                    "wb",
            ) as handle:
                pickle.dump(
                    {
                        miffi_module.GOOD: [],
                    },
                    handle,
                )

            with open(
                    logFile,
                    "w",
            ) as handle:
                handle.write("")

            class _Summary:
                def set(self, _value):
                    pass

            class _Image:
                def clone(self):
                    return self

                def getFileName(self):
                    return os.path.join(
                        tmp,
                        "mic_1.mrc",
                    )

            class _Input:
                def getSize(self):
                    return 1

                def getIdSet(self):
                    return {1}

                def __contains__(
                        self,
                        objId,
                ):
                    return objId == 1

                def getItem(
                        self,
                        field,
                        objId,
                ):
                    assert field == "id"
                    assert objId == 1
                    return _Image()

            class _Harness:
                def __init__(self):
                    self._resultsLock = (
                        threading.Lock()
                    )
                    self.processedIds = [1]
                    self.insertedIds = [1]
                    self._pendingInputIds = set()
                    self.outputCategorizeFiles = [
                        pklFile
                    ]
                    self.outputCategorizeLogFiles = [
                        logFile
                    ]
                    self.isStreamClosed = True
                    self.inputFn = "logical-input"
                    self.acceptedLabels = [
                        miffi_module.GOOD
                    ]
                    self.rejectedLabels = []
                    self.firstTime = {
                        miffi_module.OUTPUT: True,
                        miffi_module.OUTPUT_DISCARDED: True,
                    }
                    self.labelHistory = (
                        defaultdict(list)
                    )
                    self.timeHistory = []
                    self.outputLog = {}
                    self.summaryVar = _Summary()
                    self.inputSet = object()
                    self.finished = False

                def _getAllDoneIds(self):
                    return [], 0, [], []

                def _loadInputSet(
                        self,
                        _inputFn,
                ):
                    return _Input()

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass


                def _store(self):
                    pass

            protocol = _Harness()

            MiffiProtMicrographs._checkNewOutput(
                protocol
            )

            self.assertNotIn(
                1,
                protocol.processedIds,
                "An unclassified processed id must be released for retry.",
            )

            self.assertNotIn(
                1,
                protocol.insertedIds,
                "An unclassified id must become schedulable again.",
            )

            self.assertIn(
                1,
                protocol._pendingInputIds,
                "An unclassified id below the streaming watermark must return "
                "to pending input state so an open stream can retry it.",
            )

            self.assertEqual(
                [],
                protocol.outputCategorizeFiles,
                "Consumed unclassified result files must not be replayed forever.",
            )

            self.assertEqual(
                [],
                protocol.outputCategorizeLogFiles,
                "Consumed unclassified result logs must not be replayed forever.",
            )


class TestMiffiOutputIdentityRegression(unittest.TestCase):
    def testBackingFileDoesNotBecomeDurableOutputIdentity(self):
        class _FreshOutput:
            STREAM_OPEN = 1

            def __init__(self, filename=None):
                self.filename = filename
                self.loaded = False
                self.streamState = None
                self.copiedFrom = None

            def loadAllProperties(self):
                self.loaded = True
                raise AssertionError("A backing file must not restore an output that is absent from the logical protocol outputs.")

            def setStreamState(self, state):
                self.streamState = state

            def copyInfo(self, inputs):
                self.copiedFrom = inputs

        inputs = object()

        class Harness:
            def __init__(self):
                self.inputSet = _Pointer(inputs)
                self.created = 0

            def _createSetOfMicrographs(self, suffix=""):
                self.created += 1
                return _FreshOutput()

            def _getPath(self, name):
                raise AssertionError("Backing-file paths must not define MIFFI output identity.")

        protocol = Harness()
        outputSet = MiffiProtMicrographs._loadOutputSet(protocol, miffi_module.OUTPUT)

        self.assertEqual(1, protocol.created)
        self.assertFalse(outputSet.loaded)
        self.assertEqual(_FreshOutput.STREAM_OPEN, outputSet.streamState)
        self.assertIs(inputs, outputSet.copiedFrom)


class TestMiffiBatchInputIdentityRegression(unittest.TestCase):
    def testPrepareBatchPreservesMicrographsWithSameBasename(self):
        import os
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            firstDir = os.path.join(tmp, "first")
            secondDir = os.path.join(tmp, "second")
            os.makedirs(firstDir)
            os.makedirs(secondDir)

            firstMic = os.path.join(firstDir, "shared.mrc")
            secondMic = os.path.join(secondDir, "shared.mrc")

            with open(firstMic, "w") as handle:
                handle.write("first")
            with open(secondMic, "w") as handle:
                handle.write("second")

            class _Mic:
                def __init__(self, fileName):
                    self.fileName = fileName

                def clone(self):
                    return _Mic(self.fileName)

                def getFileName(self):
                    return self.fileName

            class _Input:
                def __init__(self):
                    self.items = {
                        1: _Mic(firstMic),
                        2: _Mic(secondMic),
                    }

                def __contains__(self, micId):
                    return micId in self.items

                def getItem(self, field, micId):
                    assert field == "id"
                    return self.items[micId]

            class _Harness:
                MIC_VISIBILITY_MAX_ATTEMPTS = 1
                MIC_VISIBILITY_RETRY_DELAY = 0

                def __init__(self):
                    self.input = _Input()
                    self.errors = []

                def _getTmpPath(self, name):
                    return os.path.join(tmp, name)

                def _loadInputSet(self, _unused=None):
                    return self.input

                def error(self, message):
                    self.errors.append(message)

            protocol = _Harness()

            batchDir, preparedIds = MiffiProtMicrographs._prepareBatchWithIds(
                protocol,
                [1, 2],
                1,
            )

            stagedFiles = sorted(
                name
                for name in os.listdir(batchDir)
                if name.endswith(".mrc")
            )

            self.assertEqual([1, 2], preparedIds)
            self.assertEqual(
                2,
                len(stagedFiles),
                "Distinct logical micrographs with the same source basename "
                "must not overwrite each other in a MIFFI batch workspace.",
            )


class TestMiffiBatchWorkspaceRegression(unittest.TestCase):
    def testInferenceCleansReusedBatchWorkspaceBeforeRunning(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            outDir = os.path.join(tmp, "micBatch1")
            os.makedirs(outDir)
            staleFile = os.path.join(outDir, "stale.pkl")
            with open(staleFile, "w") as handle:
                handle.write("stale")

            class Harness:
                def _getExtraPath(self, name):
                    return os.path.join(tmp, name)

                def _getInferenceParams(self, batchDir, outDir):
                    return "params"

                def runJob(self, program, params, env=None, numberOfThreads=None):
                    if os.path.exists(staleFile):
                        raise AssertionError("Reused MIFFI batch workspace must be cleaned before inference.")
                    with open(os.path.join(outDir, "fresh.pkl"), "w") as handle:
                        handle.write("fresh")

            protocol = Harness()

            with patch.object(miffi_module.Plugin, "getProgram", return_value="miffi"), patch.object(miffi_module.Plugin, "getEnviron", return_value={}):
                result = MiffiProtMicrographs.runMiffiInference(protocol, "/tmp/input-batch", 1)

            self.assertEqual(os.path.join(outDir, "fresh.pkl"), result)


    def testPrepareBatchCleansReusedTemporaryWorkspace(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            batchDir = os.path.join(tmp, "micBatch1")
            os.makedirs(batchDir)
            staleMic = os.path.join(batchDir, "stale.mrc")
            with open(staleMic, "w") as handle:
                handle.write("stale")

            inputSet = _FakeMicSet(ids=[2], tmp=tmp)
            harness = _PrepareBatchHarness(tmp, inputSet)

            preparedDir, preparedIds = MiffiProtMicrographs._prepareBatchWithIds(
                harness, [2], 1
            )

            self.assertEqual(batchDir, preparedDir)
            self.assertEqual([2], preparedIds)
            self.assertEqual(
                ["2__mic_2.mrc"],
                sorted(os.listdir(preparedDir)),
                "A reused temporary MIFFI batch directory must not retain micrographs from a previous run.",
            )


class TestMiffiCanonicalOutputRegression(unittest.TestCase):
    def testSourceRelationUsesCanonicalPublishedOutput(self):
        import pickle
        import tempfile
        import threading
        from collections import defaultdict

        with tempfile.TemporaryDirectory() as tmp:
            pklFile = os.path.join(tmp, "batch.pkl")
            logFile = os.path.join(tmp, "batch.log")
            with open(pklFile, "wb") as handle:
                pickle.dump({miffi_module.GOOD: ["mic_1.mrc"]}, handle)
            with open(logFile, "w") as handle:
                handle.write("")

            class _Summary:
                def set(self, value):
                    pass

            class _Image:
                def clone(self):
                    return self

                def getFileName(self):
                    return os.path.join(tmp, "mic_1.mrc")

            class _Input:
                def getSize(self):
                    return 1

                def getIdSet(self):
                    return {1}

                def __contains__(self, objId):
                    return objId == 1

                def getItem(self, field, objId):
                    return _Image()

            class _Output:
                def __init__(self, persistedIds=()):
                    self.persistedIds = set(persistedIds)

                def append(self, image):
                    pass

                def __contains__(self, objId):
                    return objId in self.persistedIds

            class Harness:
                def __init__(self):
                    self._resultsLock = threading.Lock()
                    self.processedIds = [1]
                    self.outputCategorizeFiles = [pklFile]
                    self.outputCategorizeLogFiles = [logFile]
                    self.isStreamClosed = True
                    self.inputFn = "logical-input"
                    self.acceptedLabels = [miffi_module.GOOD]
                    self.rejectedLabels = []
                    self.firstTime = {miffi_module.OUTPUT: True, miffi_module.OUTPUT_DISCARDED: True}
                    self.labelHistory = defaultdict(list)
                    self.timeHistory = []
                    self.outputLog = {}
                    self.summaryVar = _Summary()
                    self.inputSet = object()
                    self._inputClass = object
                    self._baseName = "micrographs.sqlite"
                    self.provisional = _Output()
                    self.canonical = _Output({1})
                    self.relationTarget = None
                    self.doneCalls = 0

                def _getAllDoneIds(self):
                    self.doneCalls += 1
                    return ([], 0, [], []) if self.doneCalls == 1 else ([1], 1, [1], [])

                def _loadInputSet(self, inputFn):
                    return _Input()

                def _loadOutputSet(self, outputName, suffix=""):
                    return self.provisional

                def _updateOutputSet(self, outputName, outputSet, streamMode):
                    setattr(self, outputName, self.canonical)

                def _defineSourceRelation(self, source, target):
                    self.relationTarget = target

                def _plotMiffiLabelHistogram(self):
                    pass

                def _plotMiffiTimeEvolution(self):
                    pass


                def _store(self):
                    pass

            protocol = Harness()

            with patch.object(miffi_module, "populate_and_update_categories", return_value={}):
                MiffiProtMicrographs._checkNewOutput(protocol)

            self.assertIs(protocol.canonical, protocol.relationTarget, "Source relations must target the canonical output published by the runtime.")


class TestMiffiScipionOutputFactoryRegression(unittest.TestCase):
    def testNewOutputsUseScipionMicrographFactoryInsteadOfBackingFileNames(self):
        class _FreshOutput:
            STREAM_OPEN = 1

            def __init__(self):
                self.streamState = None
                self.copiedFrom = None

            def setStreamState(self, state):
                self.streamState = state

            def copyInfo(self, inputs):
                self.copiedFrom = inputs

        inputs = object()

        class Harness:
            def __init__(self):
                self.inputSet = _Pointer(inputs)
                self.createdSuffixes = []

            def _createSetOfMicrographs(self, suffix=''):
                self.createdSuffixes.append(suffix)
                return _FreshOutput()

            def _getPath(self, name):
                raise AssertionError("MIFFI must not choose physical backing filenames for logical outputs.")

        protocol = Harness()

        accepted = MiffiProtMicrographs._loadOutputSet(protocol, miffi_module.OUTPUT)
        discarded = MiffiProtMicrographs._loadOutputSet(protocol, miffi_module.OUTPUT_DISCARDED, suffix="_discarded")

        self.assertEqual(["", "_discarded"], protocol.createdSuffixes)
        self.assertEqual(_FreshOutput.STREAM_OPEN, accepted.streamState)
        self.assertEqual(_FreshOutput.STREAM_OPEN, discarded.streamState)
        self.assertIs(inputs, accepted.copiedFrom)
        self.assertIs(inputs, discarded.copiedFrom)


class TestMiffiContinueRegression(unittest.TestCase):
    def testContinueSkipsPersistedIdsAndSchedulesOnlyPendingOnes(self):
        class _InputSet:
            def __init__(self):
                self.closed = False

            def getIdSet(self):
                return {1, 2, 3}

            def getSize(self):
                return 3

            def isStreamClosed(self):
                return True

            def close(self):
                self.closed = True

            def hasChangedSince(self, _lastCheck):
                return True

        class _BatchSize:
            def get(self):
                return 1

        class _Harness:
            def __init__(self):
                self._inputSet = _InputSet()
                self.inputSet = _Pointer(self._inputSet)
                self.inputFn = "logical-input"
                self.insertedIds = []
                self.isStreamClosed = False
                self.streamingBatchSize = _BatchSize()
                self.scheduled = []
                self.newDeps = []

            def debug(self, _message):
                pass

            def info(self, _message):
                pass

            def _loadInputSet(self, _inputFn):
                return self._inputSet

            def isContinued(self):
                return True

            def _getAllDoneIds(self):
                return [1, 2], 2, [1, 2], []

            def _insertNewImageSteps(self, newIds, batchSize):
                self.scheduled.append((list(newIds), batchSize))
                self.insertedIds.extend(newIds)
                return [99]

            def updateSteps(self):
                pass

        protocol = _Harness()

        MiffiProtMicrographs._checkNewInput(
            protocol
        )

        self.assertEqual(
            [([3], 1)],
            protocol.scheduled,
        )

        self.assertEqual(
            {1, 2, 3},
            set(protocol.insertedIds),
        )

        self.assertTrue(
            protocol._inputSet.closed,
        )


class TestMiffiLateVisibilityRetryRegression(unittest.TestCase):
    def testInvisibleScheduledMicReturnsToPendingAfterWatermarkAdvance(self):
        import threading

        class _Harness:
            def __init__(self):
                self.insertedIds = [7]
                self._pendingInputIds = set()
                self.processedIds = []
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self._resultsLock = threading.Lock()
                self.isStreamClosed = True
                self.deletedBatches = []

            def _prepareBatchWithIds(self, newIds, counterBatch):
                self.assertedNewIds = list(newIds)
                return "/tmp/miffi-invisible-watermarked", []

            def deleteBatch(self, batchDir):
                self.deletedBatches.append(batchDir)

            def info(self, _message):
                pass

        protocol = _Harness()

        MiffiProtMicrographs.miffStep(protocol, [7], 1)

        self.assertNotIn(
            7,
            protocol.insertedIds,
            "An invisible scheduled micrograph must be released from insertedIds.",
        )
        self.assertIn(
            7,
            protocol._pendingInputIds,
            "An invisible micrograph discovered below the advanced watermark "
            "must return to pending input state so a later poll can retry it.",
        )

    def testInvisibleMicIsNotMarkedProcessedAndCanBeScheduledAgain(self):
        import tempfile
        import threading
        from pathlib import Path

        class _InputSet:
            def __init__(self, tmp):
                self._ids = set()
                self._tmp = tmp
                self.closed = False

            def __contains__(self, micId):
                return micId in self._ids

            def getItem(self, field, micId):
                assert field == "id"
                if micId not in self._ids:
                    raise UnboundLocalError("row not found for id %r" % micId)
                return _FakeMic(micId, self._tmp)

            def getIdSet(self):
                return set(self._ids)

            def getSize(self):
                return len(self._ids)

            def isStreamClosed(self):
                return True

            def close(self):
                self.closed = True

            def hasChangedSince(self, _lastCheck):
                return True

        class _StreamingBatchSize:
            def get(self):
                return 1

        class _Harness:
            MIC_VISIBILITY_MAX_ATTEMPTS = 1
            MIC_VISIBILITY_RETRY_DELAY = 0

            def __init__(self, tmp, inputSet):
                self._tmp = tmp
                self._inputSet = inputSet
                self.inputSet = _Pointer(inputSet)
                self.inputFn = "logical-input"
                self.insertedIds = [1]
                self.processedIds = []
                self.outputCategorizeFiles = []
                self.outputCategorizeLogFiles = []
                self.isStreamClosed = True
                self._resultsLock = threading.Lock()
                self.streamingBatchSize = _StreamingBatchSize()
                self.scheduledBatches = []
                self.errors = []
                self.lastCheck = None
                self.newDeps = []

            def _getTmpPath(self, name):
                return str(Path(self._tmp) / name)

            def _loadInputSet(self, _inputFn):
                return self._inputSet

            def prepareBatch(self, newIds, counterBatch):
                return MiffiProtMicrographs.prepareBatch(self, newIds, counterBatch)

            def _prepareBatchWithIds(self, newIds, counterBatch):
                return MiffiProtMicrographs._prepareBatchWithIds(
                    self, newIds, counterBatch
                )

            def runMiffiInference(self, _batchDir, _counterBatch):
                return "inference.pkl"

            def runMiffiCategorize(self, _counterBatch, _inferenceOutputFile):
                return "categorize.pkl", "categorize.log"

            def deleteBatch(self, batchDir):
                Path(batchDir).rmdir()

            def error(self, message):
                self.errors.append(message)

            def info(self, _message):
                pass

            def debug(self, _message):
                pass

            def isContinued(self):
                return False

            def _insertNewImageSteps(self, newIds, _batchSize):
                self.scheduledBatches.append(list(newIds))
                return []

            def updateSteps(self):
                pass

        with tempfile.TemporaryDirectory() as tmp:
            inputSet = _InputSet(tmp)
            protocol = _Harness(tmp, inputSet)

            MiffiProtMicrographs.miffStep(protocol, [1], 1)

            self.assertNotIn(
                1,
                protocol.processedIds,
                "A micrograph omitted from the MIFFI batch because it is not yet visible must not be marked as processed.",
            )

            inputSet._ids = {1}
            MiffiProtMicrographs._checkNewInput(protocol)

            self.assertEqual(
                [[1]],
                protocol.scheduledBatches,
                "A temporarily invisible micrograph must become schedulable again once it is visible.",
            )


class TestMiffiStreamingArchitectureRegression(unittest.TestCase):
    def testProtocolUsesProtStreamingBaseInfrastructure(self):
        self.assertTrue(
            issubclass(MiffiProtMicrographs, ProtStreamingBase),
            "MIFFI streaming must inherit from ProtStreamingBase.",
        )
        self.assertIs(
            MiffiProtMicrographs._insertAllSteps,
            ProtStreamingBase._insertAllSteps,
            "MIFFI must use ProtStreamingBase._insertAllSteps instead of "
            "maintaining a plugin-local copy.",
        )
        self.assertIs(
            MiffiProtMicrographs.resumableStepGeneratorStep,
            ProtStreamingBase.resumableStepGeneratorStep,
            "MIFFI must use ProtStreamingBase.resumableStepGeneratorStep.",
        )
        self.assertIs(
            MiffiProtMicrographs._stepsCheck,
            ProtStreamingBase._stepsCheck,
            "MIFFI must use ProtStreamingBase._stepsCheck.",
        )


class TestMiffiStreamingGeneratorRegression(unittest.TestCase):
    def testGeneratorStopsImmediatelyWhenProtocolFails(self):
        from unittest.mock import Mock

        protocol = Mock()
        protocol.finished = False
        protocol.isFailed.return_value = True
        protocol._checkNewInput.side_effect = AssertionError(
            "A failed streaming generator must stop before polling input."
        )

        MiffiProtMicrographs.stepsGeneratorStep(protocol)

        protocol.initializeParams.assert_called_once()
        protocol._checkNewInput.assert_not_called()
        protocol._checkNewOutput.assert_not_called()
        protocol._insertFunctionStep.assert_not_called()

    def testGeneratorStopsImmediatelyWhenAlreadyFinished(self):
        # The old _stepsCheck polling callback (driven externally by the
        # executor, with a wait=True join step unlocked via STATUS_NEW) is
        # now a single resumable generator step whose while-loop condition
        # is self.finished itself.
        from unittest.mock import Mock

        protocol = Mock()
        protocol.finished = True

        MiffiProtMicrographs.stepsGeneratorStep(protocol)

        protocol.initializeParams.assert_called_once()
        protocol._checkNewInput.assert_not_called()
        protocol._checkNewOutput.assert_not_called()
        protocol._insertFunctionStep.assert_called_once()

    def testGeneratorPollsUntilFinishedThenInsertsOutputStep(self):
        from unittest.mock import Mock

        protocol = Mock()
        protocol.finished = False
        protocol.isFailed.return_value = False
        protocol._getStreamingSleepOnWait.return_value = 0

        polls = []

        def _checkNewOutput():
            polls.append(True)
            # Terminate on the second poll, like a closed input stream.
            if len(polls) == 2:
                protocol.finished = True

        protocol._checkNewOutput.side_effect = _checkNewOutput

        with patch.object(miffi_module.time, "sleep"):
            MiffiProtMicrographs.stepsGeneratorStep(protocol)

        self.assertEqual(2, protocol._checkNewInput.call_count)
        self.assertEqual(2, len(polls))
        protocol._insertFunctionStep.assert_called_once()
