# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
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


def getMiffiEnvName(version):
    return "miffi-%s" % version

V1_0_0 = "1.0.0" # MIFFI VERSION
VERSIONS = [V1_0_0]
MIFFI_DEFAULT_VER_NUM = V1_0_0
MIFFI_PIP_PACKAGE = "miffi"
DEFAULT_ENV_NAME = getMiffiEnvName(MIFFI_DEFAULT_VER_NUM)
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
MIFFI_ENV_ACTIVATION = 'MIFFI_ENV_ACTIVATION'
MIFFI_MODELS = "MIFFI_MODELS"
MIFFI_MODELS_DEFAULT = "miffi-models"
