# -*- coding: utf-8 -*-
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

"""
Describe your python module here:
This module will provide Miffi software: Cryo-EM micrograph filtering utilizing Fourier space information
"""
import os
from datetime import datetime
import copy
import time

from pyworkflow.protocol import STEPS_PARALLEL
import pyworkflow.protocol.params as params
from pyworkflow.utils import prettyTime, Message
from pyworkflow.utils.path import cleanPath, makePath, copyFile, moveFile, copyTree
from pwem.objects import Integer, SetOfMicrographs, Set
from pwem.protocols import EMProtocol, ProtPreprocessMicrographs
from pyworkflow.protocol.constants import STATUS_NEW
from pyworkflow import BETA, UPDATED, NEW, PROD

from .. import Plugin, MIFFI_MODELS

OUTPUT = 'outputMicrographs'
OUTPUT_DISCARDED = 'outputMicrographsDiscarded'

class MiffiProtMicrographs(ProtPreprocessMicrographs, EMProtocol):
    """
    Protocol to categorize micrographs
    """
    _label = 'categorize micrographs'
    _devStatus = BETA
    _possibleOutputs = {OUTPUT:SetOfMicrographs,
                        OUTPUT_DISCARDED:SetOfMicrographs}

    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        # You need a params to belong to a section:
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputSet', params.PointerParam, pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True)
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " Micassess can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParallelSection(threads=1, mpi=1)

        self._defineStreamingParams(form)
        form.getParam('streamingBatchSize').setDefault(5)
        form.getParam('streamingSleepOnWait').setDefault(5)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self.initializeParams()
        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=[], wait=True, needsGPU=False)

    def createOutputStep(self):
        self._closeOutputSet()

    def initializeParams(self):
        self.finished = False
        # Important to have both:
        self.insertedIds = []   # Contains images that have been inserted in a Step (checkNewInput).
        self.processedIds = [] # Ids to be register to output
        self.counterBatch = 1
        self.isStreamClosed = self.inputSet.get().isStreamClosed()
        # Contains images that have been processed in a Step (checkNewOutput).
        self.inputFn = self.inputSet.get().getFileName()
        self._inputClass = self.inputSet.get().getClass()
        self._inputType = self.inputSet.get().getClassName().split('SetOf')[1]
        self._baseName = '%s.sqlite' % self._inputType.lower()

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all ctfs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new images to process from the input set
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(self.inputFn))
        self.debug('Last check: %s, modification: %s'
                    % (prettyTime(self.lastCheck),
                       prettyTime(mTime)))
        # If the input.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and self.insertedIds:  # If this is empty it is dut to a static "continue" action or it is the first round
            return None

        inputSet = self._loadInputSet(self.inputFn)
        inputSetIds = inputSet.getIdSet()
        newIds = [idImage for idImage in inputSetIds if idImage not in self.insertedIds]

        self.lastCheck = datetime.now()
        self.isStreamClosed = inputSet.isStreamClosed()
        inputSet.close()

        outputStep = self._getFirstJoinStep()

        if self.isContinued() and not self.insertedIds:  # For "Continue" action and the first round
            doneIds, _, _, _ = self._getAllDoneIds()
            skipIds = list(set(newIds).intersection(set(doneIds)))
            newIds = list(set(newIds).difference(set(doneIds)))
            self.info("Skipping Images with ID: %s, seems to be done" % skipIds)
            self.insertedIds = doneIds  # During the first round of "Continue" action it has to be filled

        # Now handle the steps depending on the streaming batch size
        streamingBatchSize = self.streamingBatchSize.get()
        if len(newIds) < streamingBatchSize and not self.isStreamClosed:
            return  # No register any step if the batch size is not reach unless is the lass iter

        if newIds:
            streamingBatchSize = len(newIds) if streamingBatchSize == 0 else streamingBatchSize
            fDeps = self._insertNewImageSteps(newIds, streamingBatchSize)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        doneListIds, currentOutputSize, _, _ = self._getAllDoneIds()
        processedIds = copy.deepcopy(self.processedIds)
        newDone = [imageId for imageId in processedIds if imageId not in doneListIds]
        allDone = len(doneListIds) + len(newDone)
        maxSize = self._loadInputSet(self.inputFn).getSize()

        # We have finished when there is not more input images
        # (stream closed) or when the limit of output size is met
        self.finished = self.isStreamClosed and allDone == maxSize
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if not self.finished and not newDone:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        inputSet = self._loadInputSet(self.inputFn)
        outputSet = self._loadOutputSet(self._inputClass, self._baseName)

        for imageId in newDone:
            image = inputSet.getItem("id", imageId).clone()
            outputSet.append(image)

        self._updateOutputSet(OUTPUT, outputSet, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

        self._store()

    def _loadInputSet(self, inputFn):
        self.debug("Loading input db: %s" % inputFn)
        inputSet = self._inputClass(filename=inputFn)
        inputSet.loadAllProperties()
        return inputSet

    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)

        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputSet.get()
        outputSet.copyInfo(inputs)

        return outputSet

    def _insertNewImageSteps(self, newIds, batchSize):
        """ Insert steps to register new images (from streaming)
        Params:
            newIds: input images ids to be processed
        """
        deps = []
        # Loop through the image IDs in batches
        for i in range(0, len(newIds), batchSize):
            batchIds = newIds[i:i + batchSize]
            if len(batchIds) == batchSize or self.isStreamClosed:
                stepId = self._insertFunctionStep(self.miffStep, batchIds, self.counterBatch, needsGPU=True,
                                              prerequisites=[])
                self.counterBatch += 1
                self.insertedIds.extend(batchIds)
                deps.append(stepId)

        return deps

    def miffStep(self, newIds, counterBatch):
        """ Call miff with the appropriate parameters. """

        print('Hello miffi')
        directoryBatch = self.prepareBatch(newIds, counterBatch)
        print(directoryBatch)

        inference_output_file = self.runMiffiInference(counterBatch)
        self.runMiffiCategorize(counterBatch, inference_output_file)


        #self.appendTotalOutputStar(numPass)
        #self.copyMiffiOutput(numPass)
        self.processedIds.extend(newIds)

        if not self.isStreamClosed:
            self.delayRegister()

    def prepareBatch(self, newIds, counterBatch):
        batchDirTmp = self._getExtraPath('micBatch%d' % counterBatch)
        makePath(batchDirTmp)
        inputMicSet = self._loadInputSet(self.inputFn)
        for micId in newIds:
            mic = inputMicSet.getItem("id", micId).clone()
            micName = mic.getFileName()
            micFnOrig = os.path.abspath(micName)
            micDest = os.path.join(batchDirTmp, os.path.basename(micName))
            print(micFnOrig)
            print(micDest)
            copyFile(micFnOrig, micDest)

        return batchDirTmp

    def copyMiffiOutput(self, numPass):
        copyTree(self._getTmpPath('output%s' % numPass), self._getExtraPath('MicAssess'))

    def runMiffiInference(self, numPass):
        cwd = self._getExtraPath('micBatch%d' % numPass)
        args = ['inference '
                '--micdir %s ' % cwd,
                '-w %s ' % '*.mrc',
                '--outname %s ' % os.path.join(cwd, 'file'),
                '--g %(GPU)s ',
                '-m %s' % Plugin.getVar(MIFFI_MODELS)
                ]
        params = ' '.join(args)
        program = Plugin.getProgram('miffi')
        self.runJob(program, params, env=Plugin.getEnviron(), numberOfThreads=1)

        inference_output_file = find_file_with_ending(cwd, '.pkl')
        return inference_output_file

    def runMiffiCategorize(self, numPass, inference_output_file):
        cwd = self._getExtraPath('micBatch%d' % numPass)
        args = ['categorize '
                '-t %s ' % 'list',
                '-i %s ' % inference_output_file,
                #'--csg %s ' % os.path.join(cwd, 'csgfile'),
                '--sb ',
                '--sc' # Split all categories into high and low confidence based on a cutoff value.
                ]

        params = ' '.join(args)
        program = Plugin.getProgram('miffi')
        self.runJob(program, params, env=Plugin.getEnviron(), numberOfThreads=1)


    def deleteBatch(self):
        pass

    def delayRegister(self):
        delay = self.streamingSleepOnWait.get()
        self.info('Sleeping for delay time: %d' %delay)
        time.sleep(delay)

    # ------------------------- UTILS functions --------------------------------
    def _getAllDoneIds(self):
        doneIds = []
        acceptedIds = []
        discardedIds = []
        sizeOutput = 0

        if hasattr(self, OUTPUT):
            sizeOutput += self.outputMicrographs.getSize()
            acceptedIds.extend(list(self.outputMicrographs.getIdSet()))
            doneIds.extend(acceptedIds)

        if hasattr(self, OUTPUT_DISCARDED):
            sizeOutput += self.outputMicrographsDiscarded.getSize()
            discardedIds.extend(list(self.outputMicrographsDiscarded.getIdSet()))
            doneIds.extend(discardedIds)

        return doneIds, sizeOutput, acceptedIds, discardedIds

    def _getArgs(self, numPass):
        """ Return the list of args for the command. """
        args = ['-i %s ' % self._getExtraPath('micBatch%d' % numPass),
                '-o output%s ' % numPass,
                #'-m %s' % Plugin.getVar(CRYOASSESS_MODELS),
                '-b %d' % self.batchSize.get(),
                '--t1 %0.2f' % self.threshold.get(),
                '--t2 %0.2f' % self.threshold2.get(),
                '--threads %d' % self.numberOfThreads.get(),
                '--gpus %(GPU)s']

        return args

    def _summary(self):
        pass

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append(f"{self.message} has been printed in this run {self.times} times.")
            if self.previousCount.hasPointer():
                methods.append("Accumulated count from previous runs were %i."
                               " In total, %s messages has been printed."
                               % (self.previousCount, self.count))
        return methods


def find_file_with_ending(directory, ending):
    """
    Finds a file with the specified ending in the given directory.

    Parameters:
        directory (str): The path to the directory to search in.
        ending (str): The file ending to search for (e.g., '.txt', '.csv').

    Returns:
        str or None: The full path of the first file found with the specified ending, or None if no such file exists.
    """
    try:
        # Iterate through files in the directory
        for filename in os.listdir(directory):
            # Check if the file ends with the specified ending
            if filename.endswith(ending):
                return os.path.join(directory, filename)
        # Return None if no file with the specified ending is found
        return None
    except FileNotFoundError:
        print(f"Error: The directory '{directory}' does not exist.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None