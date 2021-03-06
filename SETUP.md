# Note about Python environment

nuskybgd is written in Python 3. You need a Python 3 installation and PyXspec
compiled against that installation. One potential source of headache is
obtaining a Python environment in which PyXspec will work. It is recommended
that you avoid using system Python (both the executable and packages), and
instead set up your own Python installation, which you have full control over,
and won't be affected by any system update.

If you have other programs that had modified your shell rc file, to modify
your PYTHONPATH environment variable, or to change the way the Python
executable is invoked (CIAO for instance), they may be incompatible and cause
problems. If you run into problems, first check that they are not the source
of the trouble by disabling them. The set up instructions here are only
guaranteed to work in a fresh installation and vanilla environment.


## Setting up a Python 3 environment using conda

The example here uses Miniconda to install Python 3 on Ubuntu 16.04 LTS, and
sets up PyXspec using that Python installation. If you already have a Conda
installation of Python 3 but have other program dependencies layered on top,
it is best to create a fresh new environment for PyXspec and nuskybgd.

https://conda.io/projects/conda/en/latest/user-guide/install/index.html

I chose Miniconda because it just has the base conda and Python and nothing
else. If you already have conda installed, go on to create a new environment.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

When prompted, install it into somewhere in your `$HOME` directory. By
default it is `$HOME/miniconda3`.

Run `conda init` once to modify your shell rc file (it will add itself to
`$PATH`). For example, if the shell is `zsh`,

```
$HOME/miniconda3/bin/conda init zsh
```

A snippet for conda's initialization will be inserted into `~/.zshrc`.

Start a new shell session for the modified shell rc to take effect.

> **Note:**
>
> The setup flow below was originally written for Python 3.5. **The minimum
> requirement is now Python 3.6.** The commands and most of the sample output
> have been updated to reflect this change.

Create an environment named `pyxspec36` with Python 3.6

```
conda create --name pyxspec36 python=3.6
```

Show list of conda environments

```
conda info --envs
```

Activate:

```
conda activate pyxspec36
```

Now there should be an additional `(pyxspec36)` in front of the command prompt.

Verify that `python` and `pip` resolve to the new installation. The `-a`
option shows others found in your `$PATH`; the first result is the one that
will be run. If you just added a new executable to the `$PATH`, it may not be
the one selected even if it is the first entry; in that case starting a new
shell session will fix it.

```
?????? $ which -a python
/home/qw/miniconda3/envs/pyxspec36/bin/python
/usr/local/bin/python
/usr/bin/python

?????? $ which -a pip
/home/qw/miniconda3/envs/pyxspec36/bin/pip
/usr/local/bin/pip

?????? $ pip list
Package        Version
-------------- ---------
certifi    2020.6.20
Django         1.10.5
ipython        7.7.0
pip        20.2.1
prompt-toolkit 2.0.9
setuptools 49.2.1.post20200807
wheel      0.34.2
```

A very fresh installation indeed? Not quite. In this user directory, some
packages (`Django`, `ipython`, `prompt-toolkit`) are actually found in
`.local/lib/python3.6/site-packages`, which is where a different Python 3.6
installation placed modules when they were installed with `pip install --user`.
This is not good (see notes below about IPython).


**Several additional checks:**

* `python3`

The nuskybgd programs are run by `#!/usr/bin/env python3`. Thus, they need to
run in an environment in which `python3` resolves to the Miniconda Python that
was just installed. Check this using `which python3`.

```
?????? $ which -a python3
/home/qw/miniconda3/envs/pyxspec36/bin/python3
/usr/bin/python3
```

An additional `python3` executable is in the `$PATH`. If that one is run
instead, more likely than not there will be problems.


* `ipython`

You may want to run some nuskybgd functions in an interactive session. In that
case, you should check that `ipython` resolves to the Miniconda Python that
was just installed.

```
?????? $ which -a ipython
/usr/local/bin/ipython
```

In this case, there was no `ipython` executable in the expected location
`/home/qw/miniconda3/envs/pyxspec36/bin`. A quick check with `pip` shows this:

```
?????? $ pip install ipython
Requirement already satisfied: ipython in ./.local/lib/python3.6/site-packages (7.x.x)
```

`pip` sees an existing `ipython` module in its search path and does not
install `ipython` into the Miniconda Python. This is not good because that
`ipython` module was placed there by some other Python's `pip` (installed with
the `--user` flag), may be incompatible, and can change at any time. It is
best to empty this `.local/lib/python3.6` directory altogether, and for
whatever needs the other installation have, perhaps migrate to a user managed
Python installation instead (so that there is no need to use the `--user`
flag).

After removing this folder, and installing `ipython` with `pip`,

```
?????? $ which -a ipython
/home/qw/miniconda3/envs/pyxspec36/bin/ipython
/usr/local/bin/ipython
```


* `jupyter`

Install Jupyter with `pip install jupyter` and check whether the correct
executable will be run.

```
?????? $ which -a jupyter
/home/qw/miniconda3/envs/pyxspec36/bin/jupyter
/usr/local/bin/jupyter
```

Although it is possible to add the Miniconda Python executable as a kernel to
another Jupyter installation, it is simplest to run the Jupyter server
installed by the Miniconda Python.


**Finally, install these packages:**

```
pip install numpy scipy astropy
```

Do this after numpy is installed:

```
pip install pyregion
```

At this point,

```
?????? $ pip list
Package            Version
------------------ -------------------
appnope            0.1.0
argon2-cffi        20.1.0
astropy            4.0.1.post1
attrs              19.3.0
backcall           0.2.0
bleach             3.1.5
certifi            2020.6.20
cffi               1.14.1
Cython             0.29.21
decorator          4.4.2
defusedxml         0.6.0
entrypoints        0.3
importlib-metadata 1.7.0
ipykernel          5.3.4
ipython            7.16.1
ipython-genutils   0.2.0
ipywidgets         7.5.1
jedi               0.17.2
Jinja2             2.11.2
jsonschema         3.2.0
jupyter            1.0.0
jupyter-client     6.1.6
jupyter-console    6.1.0
jupyter-core       4.6.3
MarkupSafe         1.1.1
mistune            0.8.4
nbconvert          5.6.1
nbformat           5.0.7
notebook           6.1.1
numpy              1.19.1
packaging          20.4
pandocfilters      1.4.2
parso              0.7.1
pexpect            4.8.0
pickleshare        0.7.5
pip                20.2.1
prometheus-client  0.8.0
prompt-toolkit     3.0.6
ptyprocess         0.6.0
pycparser          2.20
Pygments           2.6.1
pyparsing          2.4.7
pyregion           2.0
pyrsistent         0.16.0
python-dateutil    2.8.1
pyzmq              19.0.2
qtconsole          4.7.5
QtPy               1.9.0
scipy              1.5.2
Send2Trash         1.5.0
setuptools         49.2.1.post20200807
six                1.15.0
terminado          0.8.3
testpath           0.4.4
tornado            6.0.4
traitlets          4.3.3
wcwidth            0.2.5
webencodings       0.5.1
wheel              0.34.2
widgetsnbextension 3.5.1
zipp               3.1.0
```

Install additional packages with `pip` as needed.

---

## Compiling PyXspec against the Miniconda Python installation

https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/buildinstall.html#changepyvers-label

Following the instructions from HEASoft, first find the correct compilation
flags. You must check that the `python-config` (or `python3-config`) actually
resolves to the Miniconda Python. On this user account, only `python3-config`
does:

```
?????? $ which -a python-config
/usr/bin/python-config

?????? $ which -a python3-config
/home/qw/miniconda3/envs/pyxspec36/bin/python3-config
/usr/bin/python3-config
```

Accordingly, running `python3-config` will give the information needed. The
other one will not!

```
?????? $ python3-config --cflags
-I/home/qw/miniconda3/envs/pyxspec36/include/python3.6m -I/home/qw/miniconda3/envs/pyxspec36/include/python3.6m  -Wno-unused-result -Wsign-compare -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O3 -pipe  -fdebug-prefix-map==/usr/local/src/conda/- -fdebug-prefix-map==/usr/local/src/conda-prefix -fuse-linker-plugin -ffat-lto-objects -flto-partition=none -flto -DNDEBUG -fwrapv -O3 -Wall -Wstrict-prototypes

?????? $ python3-config --ldflags
-L/home/qw/miniconda3/envs/pyxspec36/lib/python3.6/config-3.6m -L/home/qw/miniconda3/envs/pyxspec36/lib -lpython3.6m -lpthread -ldl  -lutil -lrt -lm  -Xlinker -export-dynamic
```

Modify the file `heasoft-<ver>/Xspec/BUILD_DIR/hmakerc` per instructions:

```
PYTHON_INC="-I/home/qw/miniconda3/envs/pyxspec36/include/python3.6m"
PYTHON_LIB="-L/home/qw/miniconda3/envs/pyxspec36/lib/python3.6/config-3.6m -L/home/qw/miniconda3/envs/pyxspec36/lib -lpython3.6m"
```

Then follow the instructions to rebuild PyXspec.

**Note:**

You may not have `hmake` on your system. The HEASoft installation
itself has a `hmake` executable, which will be available after HEASoft has
been initialized, e.g. with `heainit`.

**Note 2:**

You may encounter an error in the final linking step of the build
(when `hmake` is run) with the following message about LTO:

```
lto1: fatal error: bytecode stream generated with LTO version 6.0 instead of the expected 4.1
```

The `g++` compiler on this system has a different link time optimizer version
vs. the one used to compile the Miniconda Python (against which PyXspec is
being compiled). In this case, the compiler on this Ubuntu 16.04 LTS is older.
The proper way to fix this would be to use a compiler from conda for HEASoft,
but that would be too much trouble. I copied the command for the linking step
just above the error message shown above, and appended ` -fno-lto` to the end
and ran it (avoids the issue by disabling link time optimization).


Check that PyXspec has been compiled successfully.

- In a new shell session, activate the Python environment

```
conda activate pyxspec36
```

- Initialize HEASoft

```
heainit
```

- Start IPython and note the message:

```
?????? $ ipython
Python 3.6.10 |Anaconda, Inc.| (default, May  7 2020, 23:06:31)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.16.1 -- An enhanced Interactive Python. Type '?' for help.

In [1]:
```

This is indeed the one that was just installed.

```
In [1]: import xspec

In [2]:
```

There is no error message.

Next, import `astropy.io.fits` and `pyregion`:

```
In [2]: import astropy.io.fits

In [3]: import pyregion

In [4]:
```

If they imported without problem, it is most likely that you are golden.

Finally, a note about the order of imports. I had not thought the order
mattered, but on my Mac, there is an issue with pyxpec and astropy using
different versions of libcfitsio. Since HEASoft is very up-to-date, PyXspec
uses the newer version of libcfitsio with more symbols. There, PyXspec must be
imported first to avoid the problem.

---
