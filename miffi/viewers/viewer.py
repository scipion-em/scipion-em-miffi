from miffi.protocols import MiffiProtMicrographs
from miffi.protocols.protocol_miffi import MIFFI_LABEL, OUTPUT, OUTPUT_DISCARDED

from pwem.viewers import ObjectView, EmProtocolViewer
from pwem.viewers.showj import *

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
import pyworkflow.protocol.params as params

import matplotlib.pyplot as plt


class MiffiViewer(EmProtocolViewer):
    """ This viewer is intended to visualize the selection made by
        the Miffi - categorize micrographs protocol.
    """
    _label = 'viewer Categorize Micrographs'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [MiffiProtMicrographs]


    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('visualizeMicrographs', params.LabelParam,
                      label="Visualize accepted Micrographs",
                      help="Visualize Micrographs with the respective Miffi label.")
        form.addParam('visualizeDiscardedMicrographs', params.LabelParam,
                      label="Visualize discarded Micrographs",
                      help="Visualize discarded Micrographs with the respective Miffi label.")
        form.addParam('visualizeLabelsHistogram', params.LabelParam,
                      label="Visualize Miffi Labels Histogram",
                      help="Visualize a histogram with miffi labels.")
        form.addParam('visualizeLabelsVsTime', params.LabelParam,
                      label="Visualize Miffi Labels Time Evolution",
                      help="Visualize accepted and discarded labels cumulative frequency over time.")

    def _getVisualizeDict(self):
        return {
                 'visualizeMicrographs': self._visualizeMicrographsF,
                 'visualizeDiscardedMicrographs' : self._visualizeDiscardedMicrographsF,
                 'visualizeLabelsHistogram': self._visualizeMiffiHistogram,
                 'visualizeLabelsVsTime': self._visualizeMiffiTimeEvolutionPlot
                }

    def _visualizeMicrographsF(self, e=None):
        return self._visualizeMicrographs(OUTPUT)

    def _visualizeDiscardedMicrographsF(self, e=None):
        return self._visualizeMicrographs(OUTPUT_DISCARDED)

    def _visualizeMiffiHistogram(self, e=None):
        self._showImage(self.protocol.getMiffiLabelHistogramPlot())

    def _visualizeMiffiTimeEvolutionPlot(self, e=None):
        self._showImage(self.protocol.getMiffiTimeEvolutionPlot())

    def _showImage(self, imagePath):
        if os.path.exists(imagePath):
            image = plt.imread(imagePath)
            height, width, _ = image.shape
            dpi = 100
            figsize = (width / dpi, height / dpi)

            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
            ax.imshow(image, aspect='auto')
            ax.axis('off')  # Hide axes completely
            for spine in ax.spines.values():
                spine.set_visible(False)  # Hide borders

            fig.patch.set_facecolor('white')  # Ensure white background
            plt.tight_layout()
            plt.show()

    def _visualizeMicrographs(self, objName):
        views = []

        labels = ('id enabled psdCorr._filename _filename plotGlobal._filename thumbnail._filename %s'
                  % MIFFI_LABEL)
        labelRender = ('psdCorr._filename plotGlobal._filename thumbnail._filename')

        if self.protocol.hasAttribute(objName):
            setMicrographs = getattr(self.protocol, objName)
            views.append(ObjectView(
                self._project, setMicrographs.getObjId(), setMicrographs.getFileName(),
                viewParams={MODE: MODE_MD,
                            ORDER: labels,
                            RENDER: labelRender,
                            VISIBLE: labels,
                            ZOOM: 5}))
        else:
            self.infoMessage('%s does not have %s%s'
                             % (self.protocol.getObjLabel(), objName,
                                getStringIfActive(self.protocol)),
                             title='Info message').show()
        return views

def getStringIfActive(prot):
    return ', yet.' if prot.isActive() else '.'