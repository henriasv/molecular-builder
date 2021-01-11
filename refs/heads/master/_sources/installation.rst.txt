Installation
========================


Prerequisites
-------------------------

- `libffi-dev` (for some reason)
- `packmol` (for packing water )

We recommend you use pyenv with `pyenv-virtualenv <https://github.com/pyenv/pyenv-virtualenv>`_, in which case you can run for instance 

.. code-block:: bash

    pyenv install 3.8.2
    pyenv virtualenv 3.8.2 molecular_builder
    pyenv activate molecular_builder

.. warning::
    Make sure you have installed `libffi` before you install python. Otherwise, you may end up with an obscure error involving `_ctypes`.
    This is done by :code:`brew install libffi` on mac or :code:`sudo apt install libffi-dev` on ubuntu. You then have to reinstall python (which is simple using pyenv).

Then you install `molecular_builder` itself.

.. code-block:: bash

    pip install molecular-builder


For installing packmol, something like this should work on unix systems (This is what we use when automatically generating the documentation on GitHub Actions):

.. code-block:: bash

    wget https://github.com/m3g/packmol/archive/20.010.zip
    unzip 20.010.zip
    cd packmol-20.010
    sed 's/\/usr\/bin\/gfortran/gfortran/g' Makefile > tmp.txt
    mv tmp.txt Makefile
    make
    sudo cp packmol /usr/local/bin/

The regex-replacement in the makefile may not be necessary, depending on your gfortran installation.
