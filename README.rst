=================
Miffi plugin
=================

This plugin provides a wrapper for `miffi <https://github.com/ando-lab/miffi?tab=readme-ov-file>`_ software tools for automatic micrograph assessment.

Miffi: Cryo-EM micrograph filtering utilizing Fourier space information



Installation
-------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-miffi

b) Developer's version

   * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-miffi.git

   * install

    .. code-block::

       scipion installp -p /path/to/scipion-em-miffi --devel

miffi software will be installed automatically with the plugin but you can also use an existing installation by providing *miffi_ENV_ACTIVATION* (see below).
You also have to download training models separately (see below).

**Important:** you need to have conda (miniconda3 or anaconda3) pre-installed to use this program.

Configuration variables
-----------------------

*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*MIFFI_ENV_ACTIVATION* (default = conda activate miffi-1.0.0):
Command to activate the miffi environment.

The deep-learning models can be downloaded from
`authors' website <https://cosmic-cryoem.org/software/cryo-assess/>`_ and the folder with models is set with:

*MIFF_MODELS* (default = software/em/miffi-models)

Verifying
---------

To check the installation, simply run the following Scipion test:

``scipion test miffi.tests.test_protocols_miffi.TestMiffi``

Supported versions
------------------

1.0.0

Protocols
----------

* categorize micrographs

References
-----------

1. High-Throughput Cryo-EM Enabled by User-Free Preprocessing Routines. Yilai Li, Jennifer N.Cash, John J.G. Tesmer, Michael A.Cianfrocco. Structure 2020, Volume 28 (7), Pages 858-869.e3
