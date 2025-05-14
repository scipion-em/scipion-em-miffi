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
import shutil
from datetime import datetime
import copy
import time
import pickle
from pathlib import Path
import re
from collections import defaultdict

from pyworkflow.protocol import STEPS_PARALLEL
import pyworkflow.protocol.params as params
from pyworkflow.utils import prettyTime, Message
from pyworkflow.utils.path import makePath, copyFile, copyTree
from pwem.objects import SetOfMicrographs, Set, String
from pwem.protocols import EMProtocol, ProtPreprocessMicrographs
from pyworkflow.protocol.constants import STATUS_NEW, LEVEL_ADVANCED
from pyworkflow import BETA, UPDATED, NEW, PROD

from .. import Plugin, MIFFI_MODELS

OUTPUT = 'outputMicrographs'
OUTPUT_DISCARDED = 'outputMicrographsDiscarded'

MIFFI_LABEL = '_miffi_label'
MIFFI_COMMENT = '_miffi_comment'
GOOD = 'good'
BAD_SINGLE = 'bad_single'
BAD_FILM = 'bad_film'
BAD_DRIFT = 'bad_drift'
BAD_MINOR_CRYSTALLINE = 'bad_minor_crystalline'
BAD_MAJOR_CRYSTALLINE = 'bad_major_crystalline'
BAD_CONTAMINATION = 'bad_contamination'
BAD_MULTIPLE = 'bad_multiple'
CATEGORIES = [GOOD, BAD_SINGLE, BAD_FILM, BAD_DRIFT, BAD_MINOR_CRYSTALLINE, BAD_MAJOR_CRYSTALLINE, BAD_CONTAMINATION, BAD_MULTIPLE]

class MiffiProtMicrographs(ProtPreprocessMicrographs, EMProtocol):
    """
    Protocol to categorize micrographs based on the image and the FT. It calls miffis inference and categorize programs.
    """
    _label = 'categorize micrographs'
    _devStatus = BETA
    _possibleOutputs = {OUTPUT:SetOfMicrographs,
                        OUTPUT_DISCARDED:SetOfMicrographs}

    LABELS = 0
    THRESHOLD = 1

    def __init__(self, **kwargs):
        ProtPreprocessMicrographs.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
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

        form.addSection(label='Output')
        form.addParam('doRejection', params.BooleanParam, default=True, label='Reject bad micrographs?',
                      help='If accepted, micrographs will be rejected based on the following miffi labels: \n'
                           'bad_film, bad_drift, bad_minor_crystalline, bad_major_crystalline, bad_contamination and bad_multiple.')
        form.addParam('removeBadFilm', params.BooleanParam, default=True, label='Reject bad film micrographs?',
                      condition='doRejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='If accepted, micrographs will be rejected based on bad_film label: \n'
                           'this label describes the degree of the support film coverage in the micrograph, \n'
                           'which can be one of the following: no film, minor film, major film and film.')
        form.addParam('removeBadDrift', params.BooleanParam, default=True, label='Reject bad drift micrographs?',
                      condition='doRejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='If accepted, micrographs will be rejected based on bad_drift label: \n'
                           'This label is binary and describes issues of sample displacement, \n'
                           'where sample drift, cracks, or an empty micrograph is labeled as bad drift.')
        form.addParam('removeBadMinorCryst', params.BooleanParam, default=True,
                      label='Reject bad minor crystalline ice micrographs?',
                      condition='doRejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='If accepted, micrographs will be rejected based on bad_minor_crystalline label: \n'
                           'This label describes the crystallinity of the ice in the sample, \n'
                           'which is determined by non-diffuse intensity at around 1/3.7 Å-1 in Fourier space. \n'
                           '_Note_: this will reject only minor crystalline ice micrographs.')
        form.addParam('removeBadMajorCryst', params.BooleanParam, default=True,
                      label='Reject bad major crystalline ice micrographs?',
                      condition='doRejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='If accepted, micrographs will be rejected based on bad_major_crystalline label: \n'
                           'This label describes the crystallinity of the ice in the sample, \n'
                           'which is determined by non-diffuse intensity at around 1/3.7 Å-1 in Fourier space. \n'
                           '_Note_: this will reject only major crystalline ice micrographs.')
        form.addParam('removeBadContamination', params.BooleanParam, default=True, label='Reject bad contaminated micrographs?',
                      condition='doRejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='If accepted, micrographs will be rejected based on bad_contamination label: \n'
                           'This label is binary and describes whether the micrograph is covered \n'
                           'largely in contaminant objects, which includes ice crystals and ethane contamination.')
        form.addParam('removeBadMultiple', params.BooleanParam, default=True,
                      label='Reject bad multiple micrographs?',
                      condition='doRejection',
                      expertLevel=LEVEL_ADVANCED,
                      help='If accepted, micrographs will be rejected based on bad_multiple label: \n'
                           'This label is binary and describes whether the micrograph is rejected \n'
                           'based of more than one criteria.')

        form.addParallelSection(threads=1)

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
        self.outputCategorizeFiles = []
        self.outputCategorizeLogFiles = []
        self.acceptedLabels, self.rejectedLabels = self._getDefinedLabels()
        self.outputLog = {}
        self.counterBatch = 1
        self.isStreamClosed = self.inputSet.get().isStreamClosed()
        # Contains images that have been processed in a Step (checkNewOutput).
        self.inputFn = self.inputSet.get().getFileName()
        self._inputClass = self.inputSet.get().getClass()
        self._inputType = self.inputSet.get().getClassName().split('SetOf')[1]
        self._baseName = '%s.sqlite' % self._inputType.lower()

    def _getDefinedLabels(self):
        categories = {
            'accept': [GOOD], # always accept 'good'
            'reject': []
        }

        if self.removeBadFilm:
            categories['reject'].append(BAD_FILM)
        else:
            categories['accept'].append(BAD_FILM)

        if self.removeBadDrift:
            categories['reject'].append(BAD_DRIFT)
        else:
            categories['accept'].append(BAD_DRIFT)

        if self.removeBadMinorCryst:
            categories['reject'].append(BAD_MINOR_CRYSTALLINE)
        else:
            categories['accept'].append(BAD_MINOR_CRYSTALLINE)

        if self.removeBadMajorCryst:
            categories['reject'].append(BAD_MAJOR_CRYSTALLINE)
        else:
            categories['accept'].append(BAD_MAJOR_CRYSTALLINE)

        if self.removeBadContamination:
            categories['reject'].append(BAD_CONTAMINATION)
        else:
            categories['accept'].append(BAD_CONTAMINATION)

        if self.removeBadMultiple:
            categories['reject'].append(BAD_MULTIPLE)
        else:
            categories['accept'].append(BAD_MULTIPLE)

        accepted_labels = categories['accept']
        rejected_labels = categories['reject']

        self.info("These are the labels that will be accepted %s" % accepted_labels)
        self.info("These are the labels that will be rejected %s" % rejected_labels)

        return accepted_labels, rejected_labels

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

        outputCategorizeFiles = copy.deepcopy(self.outputCategorizeFiles)
        outputCategorizeLogFiles = copy.deepcopy(self.outputCategorizeLogFiles)
        self.outputCategorizeFiles = []  # These mics are already registered
        self.outputCategorizeLogFiles = []
        inputSet = self._loadInputSet(self.inputFn)

        categorized_micrographs = defaultdict(list)
        accepted = {}
        rejected = {}

        for pkl_files in outputCategorizeFiles:
            with open(pkl_files, 'rb') as file:
                data = pickle.load(file)
                print(data)
                for category in CATEGORIES:
                    if category in data:
                        for path in data[category]:
                            mic_name = Path(path).name
                            categorized_micrographs[mic_name].append(category)

        for mic_name, labels in categorized_micrographs.items():
            # Handle bad_multiple separately
            if BAD_MULTIPLE in labels:
                other_reasons = [l for l in labels if l != BAD_MULTIPLE]
                if BAD_MULTIPLE in self.acceptedLabels:
                    accepted[mic_name] = {
                        'label': BAD_MULTIPLE,
                        'comment': ', '.join(other_reasons) if other_reasons else ''
                    }
                else:
                    rejected[mic_name] = {
                        'label': BAD_MULTIPLE,
                        'comment': ', '.join(other_reasons) if other_reasons else ''
                    }
            else:
                # Find which of the labels belong to accept or reject lists
                matching_accept = [l for l in labels if l in self.acceptedLabels]
                matching_reject = [l for l in labels if l in self.rejectedLabels]

                if matching_accept:
                    # Pick first accepted label, keep original meaning
                    accepted[mic_name] = {'label': matching_accept[0]}
                elif matching_reject:
                    rejected[mic_name] = {'label': matching_reject[0]}

        # Load output sets
        if accepted:
            outputSet = self._loadOutputSet(self._inputClass, self._baseName)
        if rejected:
            outputSetDiscarded = self._loadOutputSet(self._inputClass, 'micrographDISCARDED.sqlite')

        # Assign micrographs to their sets with attributes
        for imageId in newDone:
            image = inputSet.getItem("id", imageId).clone()
            micName = os.path.basename(image.getFileName())

            if micName in accepted:
                setLabel(image, MIFFI_LABEL, accepted[micName]['label'])
                if 'comment' in accepted[micName]:
                    setLabel(image, MIFFI_COMMENT, accepted[micName]['comment'])
                else:
                    setLabel(image, MIFFI_COMMENT, '')
                outputSet.append(image)

            elif micName in rejected:
                setLabel(image, MIFFI_LABEL, rejected[micName]['label'])
                if 'comment' in rejected[micName]:
                    setLabel(image, MIFFI_COMMENT, rejected[micName].get('comment', ''))
                else:
                    setLabel(image, MIFFI_COMMENT, '')
                outputSetDiscarded.append(image)

        if accepted:
            self._updateOutputSet(OUTPUT, outputSet, streamMode)
        if rejected:
            self._updateOutputSet(OUTPUT_DISCARDED, outputSetDiscarded, streamMode)

        # Summary log and display the results
        outputLogTmp = populate_and_update_categories(self.outputLog, outputCategorizeLogFiles)
        dict_str = '\n'.join(f'{k}: {v}' for k, v in outputLogTmp.items())
        self.summaryVar.set(dict_str)
        self.outputLog  = outputLogTmp

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
        batchDirTmp = self.prepareBatch(newIds, counterBatch)
        try:
            inference_output_file = self.runMiffiInference(batchDirTmp, counterBatch)
            categorize_output_file, categorize_output_log_file  = self.runMiffiCategorize(counterBatch, inference_output_file)
            self.outputCategorizeFiles.append(categorize_output_file)
            self.outputCategorizeLogFiles.append(categorize_output_log_file)
        except Exception as e:
            self.info('Batch number %d had problems with miffi execution' % counterBatch)
            self.info(e)

        self.processedIds.extend(newIds)
        # To have a control in the size of the protocol
        self.deleteBatch(batchDirTmp)

        if not self.isStreamClosed:
            self.delayRegister()

    def prepareBatch(self, newIds, counterBatch):
        batchDirTmp = self._getTmpPath('micBatch%d' % counterBatch)
        makePath(batchDirTmp)
        inputMicSet = self._loadInputSet(self.inputFn)
        for micId in newIds:
            mic = inputMicSet.getItem("id", micId).clone()
            micName = mic.getFileName()
            micFnOrig = os.path.abspath(micName)
            micDest = os.path.join(batchDirTmp, os.path.basename(micName))
            copyFile(micFnOrig, micDest)

        return batchDirTmp

    def copyMiffiOutput(self, numPass):
        copyTree(self._getTmpPath('output%s' % numPass), self._getExtraPath('MicAssess'))

    def runMiffiInference(self, batchDir, numPass):
        outDir = self._getExtraPath('micBatch%d' % numPass)
        makePath(outDir)
        params = self._getInferenceParams(batchDir, outDir)
        program = Plugin.getProgram('miffi')
        self.runJob(program, params, env=Plugin.getEnviron(), numberOfThreads=1)
        inference_output_file = find_file_with_ending(outDir, '.pkl')

        return inference_output_file

    def runMiffiCategorize(self, numPass, inference_file):
        cwd = self._getExtraPath('micBatch%d' % numPass)
        params = self._getCategorizeParams(inference_file)
        program = Plugin.getProgram('miffi')
        self.runJob(program, params, env=Plugin.getEnviron(), numberOfThreads=1)

        categorize_output_file = find_file_with_ending(cwd, 'dict.pkl')
        categorize_output_log_file = find_file_with_ending(cwd, '.log')

        return categorize_output_file, categorize_output_log_file

    def deleteBatch(self, tmpBatchDir):
        shutil.rmtree(tmpBatchDir)

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

    def _getInferenceParams(self, batchDir, outDir):
        """ Return the list of args for the command. """
        args = ['inference '
                '--micdir %s ' % batchDir,
                '-w %s ' % '*.mrc',
                '--outname %s ' % os.path.join(outDir, 'file'),
                '--g %(GPU)s ',
                '-m %s' % Plugin.getVar(MIFFI_MODELS)
                ]
        params = ' '.join(args)
        return params

    def _getCategorizeParams(self, inference_file):
        """ Return the list of args for the command. """
        args = ['categorize '
                '-t %s ' % 'list',
                '-i %s ' % inference_file,
                # '--csg %s ' % os.path.join(cwd, 'csgfile'),
                '--sb ', # Split all singly bad micrographs into individual categories. If this argument is not used, only film and minor crystalline categories will be written.
                #'--sc'  # Split all categories into high and low confidence based on a cutoff value.
                ]
        params = ' '.join(args)
        return params

    def _summary(self):
        summary = []
        summary.append(self.summaryVar.get())
        return summary

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

def populate_and_update_categories(existing_totals, log_files):
    for log_file_path in log_files:
        with open(log_file_path, 'r') as file:
            for line in file:
                # Updated regex to match category name without dots
                match = re.search(r'\)\s*([\w_]+)\.*\s*(\d+)', line)
                if match:
                    category = match.group(1)
                    count = int(match.group(2))
                    # Add new categories or update existing totals
                    existing_totals[category] = existing_totals.get(category, 0) + count
    return existing_totals

def setLabel(obj, label, value):
    if value is None:
        return
    setattr(obj, label, getScipionObj(value))

def getScipionObj(value):
    if isinstance(value, str):
        return String(value)
    else:
        return None