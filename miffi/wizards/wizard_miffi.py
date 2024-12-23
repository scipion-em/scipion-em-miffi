from pyworkflow.wizard import Wizard

from miffi.protocols import MiffiProtMicrographs


class MiffiProtMicsWizard(Wizard):
    # Dictionary to target protocol parameters
    _targets = [(MiffiProtMicrographs, ['message'])]
