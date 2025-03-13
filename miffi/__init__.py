# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# * your institution
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

import pyworkflow.utils as pwutils
import pwem
import os

from miffi.constants import *

__version__ = "1.0.0"  # plugin version
_logo = "miffi_logo.png"
_references = ['DaXu2024']


class Plugin(pwem.Plugin):
    _pathVars = [MIFFI_MODELS]
    _url = "https://github.com/scipion-em/scipion-em-miffi"
    _supportedVersions = VERSIONS  # binary version
    _homeVar = MIFFI_HOME

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(MIFFI_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)
        cls._defineEmVar(MIFFI_HOME, f"miffi-{MIFFI_DEFAULT_VER_NUM}")
        models_home = os.path.join(cls.getVar(MIFFI_HOME), MIFFI_MODELS)
        cls._defineVar(MIFFI_MODELS, models_home)

    @classmethod
    def getMiffiEnvActivation(cls):
        return cls.getVar(MIFFI_ENV_ACTIVATION)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch miffi. """
        environ = pwutils.Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addMiffiPackage(env, ver,
                                     default= ver == MIFFI_DEFAULT_VER_NUM)

    @classmethod
    def addMiffiPackage(cls, env, version, default=False):
        MIFFI_INSTALLED = '%s_%s_installed' % (MIFFI_PIP_PACKAGE, MIFFI_DEFAULT_VER_NUM)
        ENV_NAME = getMiffiEnvName(version)

        # try to get CONDA activation command
        installCmd = [cls.getCondaActivationCmd()]

        # Create the environment
        installCmd.append('conda create -y -n %s python=3.9' % ENV_NAME)

        # Activate the new environment
        installCmd.append('&& conda activate %s' % ENV_NAME)

        installCmd.append('&& pip install %s' % MIFFI_PIP_PACKAGE)

        models_home = cls.getVar(MIFFI_MODELS)
        pwutils.makePath(models_home)
        # install models
        installCmd.append('&& %s download -d %s' % (MIFFI_PIP_PACKAGE, models_home))

        # Flag installation finished
        installCmd.append('&& touch %s' % MIFFI_INSTALLED)

        miffi_commands = [(" ".join(installCmd), MIFFI_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None

        env.addPackage(MIFFI_PIP_PACKAGE, version=MIFFI_DEFAULT_VER_NUM,
                       tar='void.tgz',
                       commands=miffi_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def getProgram(cls, program):
        """ Create miffi command line. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(),
                                       cls.getMiffiEnvActivation(),
                                       program)
        return fullProgram

    # @classmethod
    # def runFidder(cls, protocol, args, cwd=None, numberOfMpi=1):
    #     """ Run fidder command from a given protocol. """
    #     cmd = cls.getCondaActivationCmd() + " "
    #     cmd += cls.getFidderEnvActivation()
    #     cmd += f" && CUDA_VISIBLE_DEVICES=%(GPU)s {FIDDER} "
    #     protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd, numberOfMpi=numberOfMpi)